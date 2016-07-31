#include <cstdio>
#include <cassert>
#include "base.h"
#include "parsing.h"
#include "structure.h"
#include "fold.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  assert(argc == 2);
  Init();
  rna_t rna = parsing::ParseRnaFromString(argv[1]);
  energy_t e = fold::Fold(rna);
  printf("Energy: %d\n%s\n", e, parsing::DotBracketFromPairs(p).c_str());
}
