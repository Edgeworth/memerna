#include <cstdio>
#include <random>
#include <chrono>
#include "bridge/bridge.h"
#include "parsing.h"

using namespace memerna;

rna_t GenerateRandomRna(int length) {
  rna_t rna(std::size_t(length), 0);
  for (int i = 0; i < length; ++i)
    rna[i] = rand() % 4;
  return rna;
}

void FuzzRna(const rna_t& rna, bool use_random_energy_model,
    const std::vector<bridge::Memerna>& memernas, const bridge::Rnastructure& rnastructure) {
  // Initialise everything here because we use goto later.
  int seed = rand();
  bool use_brute = rna.size() <= 22;
  bool dp_table_diff = false;
  int N = int(rna.size());
  dp_state_t rnastructure_state;
  folded_rna_t rnastructure_frna, brute_frna;
  energy_t rnastructure_efn = 0, brute_efn = 0;
  int st = 0, en = 0, a = 0;

  if (use_random_energy_model)
    LoadRandomEnergyModel(seed);
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

  for (int i = 0; i < int(memernas.size()); ++i) {
    if (memerna_folds[0].energy != memerna_folds[i].energy ||
        memerna_folds[0].energy != memerna_efns[i])
      goto print_diff;
  }

  if (use_brute) {
    brute_frna = fold::FoldBruteForce(rna, nullptr);
    if (memerna_folds[0].energy != brute_frna.energy)
      goto print_diff;
  }

  if (!use_random_energy_model) {
    rnastructure_frna = rnastructure.FoldAndDpTable(rna, &rnastructure_state);
    rnastructure_efn = rnastructure.Efn(rnastructure_frna);
    if (memerna_folds[0].energy != rnastructure_frna.energy ||
        memerna_folds[0].energy != rnastructure_efn)
      goto print_diff;
  }

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
            goto print_diff;
          }
        }
        energy_t rnastructureval = constants::MAX_E;
        if (a == fold::DP_P)
          rnastructureval = rnastructure_state.v.f(st + 1, en + 1);
        else if (a == fold::DP_U)
          rnastructureval = rnastructure_state.w.f(st + 1, en + 1);
        if (rnastructureval != constants::MAX_E && !use_random_energy_model) {
          if ((memerna0 < constants::CAP_E) != (rnastructureval < INFINITE_ENERGY - 1000) ||
              (memerna0 < constants::CAP_E && memerna0 != rnastructureval)) {
            dp_table_diff = true;
            goto print_diff;
          }
        }
      }
    }
  }
  return;
  print_diff:
  printf("Difference on len %d RNA %s\n", int(rna.size()), parsing::RnaToString(rna).c_str());
  if (use_random_energy_model)
    printf("  Using random energy model with seed: %d\n", seed);
  else
    printf("  Using T04 energy model.\n");
  for (int i = 0; i < int(memernas.size()); ++i) {
    printf("  Fold%d: %d (dp), %d (efn) - %s\n", i, memerna_folds[i].energy, memerna_efns[i],
        parsing::PairsToDotBracket(memerna_folds[i].p).c_str());
  }
  if (use_brute)
    printf("  BruteFold: %d (dp), %d (efn) - %s\n", brute_frna.energy, brute_efn,
        parsing::PairsToDotBracket(brute_frna.p).c_str());
  if (!use_random_energy_model)
    printf("  RNAstructure: %d (dp), %d (efn) - %s\n", rnastructure_frna.energy,
        rnastructure_efn, parsing::PairsToDotBracket(rnastructure_frna.p).c_str());
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

int main(int argc, char* argv[]) {
  srand(static_cast<unsigned int>(time(NULL)));
  ArgParse argparse({
      {"print-interval", ArgParse::option_t("status update every n seconds").Arg("-1")},
      {"random", ArgParse::option_t("use random energy models (disables comparison to RNAstructure)")}
  });
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(
      pos.size() == 2,
      "require size and variance");
  int base_len = atoi(pos[0].c_str());
  int variance = atoi(pos[1].c_str());
  verify_expr(base_len > 0, "invalid length");
  verify_expr(variance >= 0, "invalid variance");

  bridge::Rnastructure rnastructure("extern/rnark/data_tables/", false);
  std::vector<bridge::Memerna> memernas;
  for (const auto& fold_fn : fold::FOLD_FUNCTIONS) {
    memernas.emplace_back("data/", fold_fn);
  }

  auto start_time = std::chrono::steady_clock::now();
  auto interval = atoi(argparse.GetOption("print-interval").c_str());
  for (int i = 0;; ++i) {
    int length = base_len;
    if (variance) length += rand() % variance;
    if (interval > 0 && std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - start_time).count() > interval) {
      printf("Fuzzed %d RNA\n", i);
      start_time = std::chrono::steady_clock::now();
    }
    auto rna = GenerateRandomRna(length);
    FuzzRna(rna, argparse.HasFlag("random"), memernas, rnastructure);
  }
}

