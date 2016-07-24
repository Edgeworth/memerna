#include <cstdio>
#include <cassert>
#include "base.h"
#include "parsing.h"
#include "structure.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  assert(argc == 3);
  Init();
  folded_rna_t frna = parsing::ParseDotBracketRna(argv[1], argv[2]);
  printf("Computing energy.\n");
  std::unique_ptr<structure::Structure> structure;
  printf("Computed energy: %lf\n", energy::ComputeEnergy(frna, &structure) / 10.0);
  auto descs = structure->Description();
  for (const auto& desc : descs) {
    printf("%s\n", desc.c_str());
  }
}
