#include <cstdio>
#include <cassert>
#include <random>
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
        "Diff. DP: %d !=  EFN: %d, on: %s\n",
        dp, efn, parsing::StringFromRna(rna).c_str());
    return;
  }
}

int main() {
  Init();
  srand(time(NULL));
//  for (int i = 1; i <= 15; ++i) {
//    printf("Fuzzing len %d\n", i);
//    FuzzRnaOfLength(i);
//  }
  for (int i = 0; i < 100000; ++i) {
    int length = 200 + rand() % 20;
    if (i % 10000 == 0) printf("Fuzzed %d RNA\n", i);
    FuzzRandomRna(length);
  }
}
