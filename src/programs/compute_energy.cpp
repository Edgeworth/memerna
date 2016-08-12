#include <cstdio>
#include <cassert>
#include "base.h"
#include "parsing.h"
#include "structure.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  verify_expr(argc == 3, "requires two arguments");
  LoadEnergyModelFromDataDir();
  auto frna = parsing::ParseDotBracketRna(argv[1], argv[2]);
  std::unique_ptr<structure::Structure> structure;
  printf("Energy: %d\n", energy::ComputeEnergy(frna, &structure));
  auto descs = structure->Description();
  for (const auto& desc : descs) {
    printf("%s\n", desc.c_str());
  }
}
