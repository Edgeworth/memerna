#include <cstdio>
#include <random>
#include <chrono>
#include "bridge/bridge.h"
#include "parsing.h"
#include "energy/energy_model.h"

using namespace memerna;

template<typename RandomEngine>
void FuzzRna(const rna_t& rna, bool use_random_energy_model,
    const std::vector<bridge::Memerna>& memernas, const bridge::Rnastructure& rnastructure,
    RandomEngine& eng) {
  uint32_t seed = eng();
  if (use_random_energy_model)
    energy::LoadRandomEnergyModel(seed);

  std::vector<fold::fold_state_t> memerna_states;
  std::vector<folded_rna_t> memerna_folds;
  std::vector<energy_t> memerna_efns;
  for (const auto& memerna : memernas) {
    fold::fold_state_t state;
    auto frna = memerna.FoldAndDpTable(rna, &state);
    memerna_efns.push_back(energy::ComputeEnergy(frna));
    memerna_states.push_back(std::move(state));
    memerna_folds.push_back(std::move(frna));
  }

  bool mfe_diff = false;
  for (int i = 0; i < int(memernas.size()); ++i) {
    if (memerna_folds[0].energy != memerna_folds[i].energy ||
        memerna_folds[0].energy != memerna_efns[i])
      mfe_diff = true;
  }

  bool use_brute = rna.size() <= 24;
  folded_rna_t brute_frna;
  if (use_brute) {
    brute_frna = fold::FoldBruteForce(rna, nullptr);
    if (memerna_folds[0].energy != brute_frna.energy)
      mfe_diff = true;
  }

  folded_rna_t rnastructure_frna;
  energy_t rnastructure_efn = 0;
  dp_state_t rnastructure_state;
  if (!use_random_energy_model) {
    rnastructure_frna = rnastructure.FoldAndDpTable(rna, &rnastructure_state);
    rnastructure_efn = rnastructure.Efn(rnastructure_frna);
    if (memerna_folds[0].energy != rnastructure_frna.energy ||
        memerna_folds[0].energy != rnastructure_efn)
      mfe_diff = true;
  }


  int st = 0, en = 0, a = 0;
  bool dp_table_diff = false;
  int N = int(rna.size());
  for (st = N - 1; st >= 0; --st) {
    for (en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      for (a = 0; a < fold::DP_SIZE; ++a) {
        auto memerna0 = memerna_states[0].dp_table[st][en][a];
        for (int i = 0; i < int(memernas.size()); ++i) {
          auto memernai = memerna_states[i].dp_table[st][en][a];
          // If meant to be infinity and not.
          if (((memerna0 < constants::CAP_E) != (memernai < constants::CAP_E)) ||
              (memerna0 < constants::CAP_E && memerna0 != memernai)) {
            dp_table_diff = true;
            goto loop_end;
          }
        }
        if (!use_random_energy_model) {
          energy_t rnastructureval = constants::MAX_E;
          if (a == fold::DP_P)
            rnastructureval = rnastructure_state.v.f(st + 1, en + 1);
          else if (a == fold::DP_U)
            rnastructureval = rnastructure_state.w.f(st + 1, en + 1);
          if (rnastructureval != constants::MAX_E &&
              ((memerna0 < constants::CAP_E) != (rnastructureval < INFINITE_ENERGY - 1000) ||
              (memerna0 < constants::CAP_E && memerna0 != rnastructureval))) {
            dp_table_diff = true;
            goto loop_end;
          }
        }
      }
    }
  }
  loop_end:;
  if (mfe_diff || dp_table_diff) {
    printf("Difference on len %d RNA %s\n", int(rna.size()), parsing::RnaToString(rna).c_str());
    if (use_random_energy_model)
      printf("  Using random energy model with seed: %d\n", seed);
    else
      printf("  Using T04 energy model.\n");
    if (mfe_diff) {
      for (int i = 0; i < int(memernas.size()); ++i) {
        printf("  Fold%d: %d (dp), %d (efn) - %s\n", i, memerna_folds[i].energy, memerna_efns[i],
            parsing::PairsToDotBracket(memerna_folds[i].p).c_str());
      }
      if (use_brute)
        printf("  BruteFold: %d - %s\n", brute_frna.energy, parsing::PairsToDotBracket(brute_frna.p).c_str());
      if (!use_random_energy_model)
        printf("  RNAstructure: %d (dp), %d (efn) - %s\n", rnastructure_frna.energy,
            rnastructure_efn, parsing::PairsToDotBracket(rnastructure_frna.p).c_str());
    }
    if (dp_table_diff) {
      printf("  DP table difference at %d %d %d:\n", st, en, a);
      for (int i = 0; i < int(memernas.size()); ++i)
        printf("    Fold%d: %d\n", i, memerna_states[i].dp_table[st][en][a]);
      if (!use_random_energy_model)
        printf("    RNAstructure: V: %d W: %d\n",
            rnastructure_state.v.f(st + 1, en + 1),
            rnastructure_state.w.f(st + 1, en + 1));
    }
  }
}

int main(int argc, char* argv[]) {
  std::mt19937 eng(uint32_t(time(nullptr)));
  ArgParse argparse({
      {"print-interval", ArgParse::option_t("status update every n seconds").Arg("-1")},
      {"random", ArgParse::option_t("use random energy models (disables comparison to RNAstructure)")}
  });
  argparse.AddOptions(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(
      pos.size() == 2,
      "require min and max length");
  int min_len = atoi(pos[0].c_str());
  int max_len = atoi(pos[1].c_str());
  verify_expr(min_len > 0, "invalid min length");
  verify_expr(max_len >= min_len, "invalid max len");
  energy::LoadEnergyModelFromArgParse(argparse);

  bridge::Rnastructure rnastructure("extern/rnark/data_tables/", false);
  std::vector<bridge::Memerna> memernas;
  for (const auto& fold_fn : fold::FOLD_FUNCTIONS) {
    memernas.emplace_back(fold_fn);
  }

  auto start_time = std::chrono::steady_clock::now();
  auto interval = atoi(argparse.GetOption("print-interval").c_str());
  std::uniform_int_distribution<int> len_dist(min_len, max_len);
  for (int i = 0;; ++i) {
    int length = len_dist(eng);
    if (interval > 0 && std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - start_time).count() > interval) {
      printf("Fuzzed %d RNA\n", i);
      start_time = std::chrono::steady_clock::now();
    }
    auto rna = GenerateRandomRna(length, eng);
    FuzzRna(rna, argparse.HasFlag("random"), memernas, rnastructure, eng);
  }
}

