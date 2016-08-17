#include <cstdio>
#include <cassert>
#include <random>
#include <chrono>
#include <bridge/bridge.h>
#include "argparse.h"
#include "parsing.h"
#include "fold/fold.h"

using namespace memerna;

void FuzzRnaOfLength(int length) {
  assert(length * 2 < 31);
  rna_t rna(std::size_t(length), 0);
  for (int i = 0; i < 1 << (length * 2); ++i) {
    int data = i;
    for (int j = 0; j < length; ++j) {
      rna[j] = base_t(data & 0b11);
      data >>= 2;
    }
    SetRna(rna);
    energy_t a = fold::Fold();
    energy_t a_efn = energy::ComputeEnergy();
    std::string a_db = parsing::PairsToDotBracket(p).c_str();
    energy_t b = fold::FoldBruteForce();
    energy_t b_efn = energy::ComputeEnergy();
    std::string b_db = parsing::PairsToDotBracket(p).c_str();

    if (a != b || a != a_efn || b != b_efn) {
      printf(
          "Diff. DP: %d (efn %d) != BF: %d (efn %d), on: %s\n%s vs %s\n",
          a, a_efn, b, b_efn, parsing::RnaToString(rna).c_str(), a_db.c_str(), b_db.c_str());
      return;
    }
  }
}

rna_t GenerateRandomRna(int length) {
  rna_t rna(std::size_t(length), 0);
  for (int i = 0; i < length; ++i)
    rna[i] = rand() % 4;
  return rna;
}

void FuzzRna(const rna_t& rna,
    const bridge::Memerna& memerna, const bridge::Rnastructure& rnastructure) {
  auto memerna_dp = memerna.Fold(rna);
  auto memerna_efn = memerna.Efn(memerna_dp);
  auto rnastructure_dp = rnastructure.Fold(rna);
  auto rnastructure_efn = rnastructure.Efn(rnastructure_dp);

  if (memerna_dp.energy != memerna_efn ||
      memerna_dp.energy != rnastructure_dp.energy ||
      memerna_dp.energy != rnastructure_efn) {
    printf(
        "Diff on %s\n  Rnastructure: dp %d, efn %d\n    %s\n  Memerna: dp %d, efn %d\n    %s\n\n",
        parsing::RnaToString(rna).c_str(), rnastructure_dp.energy, rnastructure_efn,
        parsing::PairsToDotBracket(rnastructure_dp.p).c_str(),
        memerna_dp.energy, memerna_efn, parsing::PairsToDotBracket(memerna_dp.p).c_str()
    );
    return;
  }
}

#include "fold/slow_fold.h"
#include "fold/fold1.h"
#include "fold/fold2.h"

void FuzzComputeTables(const rna_t& rna) {
  SetRna(rna);
  auto table1 = fold::ComputeTablesSlow();
  auto table2 = fold::ComputeTables2();
  int N = int(rna.size());
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      for (int a = 0; a < fold::DP_SIZE; ++a) {
        if (table1[st][en][a] != table2[st][en][a] && table1[st][en][a] < constants::CAP_E) {
          printf("Diff on %s\n %d %d %d: %d (good) != %d\n",
              parsing::RnaToString(rna).c_str(), st, en, a, table1[st][en][a], table2[st][en][a]);
          return;
        }
      }
    }
  }
}

// test types:
// against efn / brute fold, against rnastructure
// if against itself, try random energy models
// input types: stdin, exhaustive rnas (up to size x), random rnas (of size, and variance)
// usage: -e?, -f --stdin --brute --rand
// if brute force, up to size x
// if random, size, variance
// for efn, try all foldings as well?

int main(int argc, char* argv[]) {
  ArgParse argparse({
      {"print-interval", ArgParse::option_t("status update every n seconds").Arg("60")}
  });
  auto ret = argparse.Parse(argc, argv);
  verify_expr(
      ret.size() == 0,
      "%s\n%s\n", ret.c_str(), argparse.Usage().c_str());

  LoadEnergyModelFromDataDir("data");
  srand(static_cast<unsigned int>(time(NULL)));
//  if (choice == "1") {
//    for (int i = 1; i <= 15; ++i) {
//      printf("Fuzzing len %d\n", i);
//      FuzzRnaOfLength(i);
//    }
//  } else if (choice == "2") {
  auto pos = argparse.GetPositional();
  verify_expr(
      argparse.GetPositional().size() == 2,
      "require size and variance");
  int base_len = atoi(pos[0].c_str());
  int variance = atoi(pos[1].c_str());
  verify_expr(base_len > 0, "invalid length");
  verify_expr(variance >= 0, "invalid variance");

  bridge::Memerna memerna("data/");
  bridge::Rnastructure rnastructure("extern/rnark/data_tables/", false);

  auto start_time = std::chrono::steady_clock::now();
  auto interval = atoi(argparse.GetOption("print-interval").c_str());
  for (int i = 0;; ++i) {
    int length = base_len;
    if (variance) length += rand() % variance;
    if (std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - start_time).count() > interval) {
      printf("Fuzzed %d RNA\n", i);
      start_time = std::chrono::steady_clock::now();
    }
    auto rna = GenerateRandomRna(length);
//    FuzzRna(length, memerna, rnastructure);
    FuzzComputeTables(rna);
  }
}

