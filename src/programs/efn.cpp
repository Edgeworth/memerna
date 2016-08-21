#include <cstdio>
#include "bridge/bridge.h"
#include "parsing.h"
#include "energy/structure.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse;
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 2, "requires primary sequence and dot bracket");

  LoadEnergyModelFromDataDir("data/");
  auto frna = parsing::ParseDotBracketRna(pos.front(), pos.back());
  std::unique_ptr<structure::Structure> structure;
  printf("Energy: %d\n", energy::ComputeEnergy(frna, &structure));
  auto descs = structure->Description();
  for (const auto& desc : descs) {
    printf("%s\n", desc.c_str());
  }
}
