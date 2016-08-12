#include <cstdio>
#include <cassert>
#include <random>
#include <chrono>
#include "base.h"
#include "parsing.h"
#include "fold.h"

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
    std::string a_db = parsing::DotBracketFromPairs(p).c_str();
    energy_t b = fold::FoldBruteForce();
    energy_t b_efn = energy::ComputeEnergy();
    std::string b_db = parsing::DotBracketFromPairs(p).c_str();

    if (a != b || a != a_efn || b != b_efn) {
      printf(
          "Diff. DP: %d (efn %d) != BF: %d (efn %d), on: %s\n%s vs %s\n",
          a, a_efn, b, b_efn, parsing::StringFromRna(rna).c_str(), a_db.c_str(), b_db.c_str());
      return;
    }
  }
}

void FuzzRandomRna(int length) {
  rna_t rna(std::size_t(length), 0);
  for (int i = 0; i < length; ++i)
    rna[i] = rand() % 4;
  SetRna(rna);
  energy_t dp = fold::Fold();
  energy_t efn = energy::ComputeEnergy();

  if (dp != efn) {
    printf(
        "Diff. DP: %d !=  EFN: %d, on:\n  %s\n  %s\n",
        dp, efn, parsing::StringFromRna(rna).c_str(),
        parsing::DotBracketFromPairs(p).c_str());
    return;
  }
}

int main(int argc, char* argv[]) {
  LoadEnergyModelFromDataDir();
  srand(time(NULL));
  verify_expr(argc >= 2, "require selection; 1 == brute force, 2 == random rnas");
  std::string choice = argv[1];
  if (choice == "1") {
    for (int i = 1; i <= 15; ++i) {
      printf("Fuzzing len %d\n", i);
      FuzzRnaOfLength(i);
    }
  } else if (choice == "2") {
    int base_len = 100;
    int variance = 20;
    if (argc >= 3)
      base_len = atoi(argv[2]);
    if (argc >= 4)
      variance = atoi(argv[3]);
    verify_expr(base_len > 0, "invalid length");
    verify_expr(variance >= 0, "invalid variance");

    auto start_time = std::chrono::steady_clock::now();
    for (int i = 0; ; ++i) {
      int length = base_len;
      if (variance) length += rand() % variance;
      if (std::chrono::duration_cast<std::chrono::seconds>(
          std::chrono::steady_clock::now() - start_time).count() > 10.0) {
        printf("Fuzzed %d RNA\n", i);
        start_time = std::chrono::steady_clock::now();
      }
      FuzzRandomRna(length);
    }
  } else {
    verify_expr(false, "unknown choice");
  }
}
