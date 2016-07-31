#include <cstdio>
#include <cassert>
#include "base.h"
#include "parsing.h"
#include "fold.h"

using namespace memerna;

void FuzzRna(int length) {
  assert(length * 2 < 31);
  rna_t rna(std::size_t(length), 0);
  for (int i = 0; i < 1 << (length * 2); ++i) {
    int data = i;
    for (int j = 0; j < length; ++j) {
      rna[j] = base_t(data & 0b11);
      data >>= 2;
    }
    energy_t a = fold::Fold(rna);
    energy_t b = fold::FoldBruteForce(rna);
    if (a != b) {
      printf(
          "Error %d != %d, on:%s\n%s\n",
          a, b, parsing::StringFromRna(rna).c_str(),
          parsing::DotBracketFromPairs(p).c_str());
      return;
    }
  }
}

int main() {
  Init();
  for (int i = 1; i <= 15; ++i) {
    printf("Fuzzing len %d\n", i);
    FuzzRna(i);
  }
}
