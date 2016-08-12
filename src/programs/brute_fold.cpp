#include <cstdio>
#include <cassert>
#include "base.h"
#include "parsing.h"
#include "structure.h"
#include "fold.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  verify_expr(argc == 2, "requires 1 argument");
  LoadEnergyModelFromDataDir();
  std::unique_ptr<structure::Structure> structure;
  auto rna = parsing::ParseRnaFromString(argv[1]);
  printf("Energy: %d\n", fold::FoldBruteForce(rna, &structure));
  printf("DB: %s\n", parsing::DotBracketFromPairs(p).c_str());
  auto descs = structure->Description();
  for (const auto& desc : descs) {
    printf("%s\n", desc.c_str());
  }
}
