#include <cstdio>
#include "parsing.h"
#include "energy/structure.h"
#include "energy/load_model.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  const auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 2, "requires primary sequence and dot bracket");

  const auto em = energy::LoadEnergyModelFromArgParse(argparse);
  const auto secondary = parsing::ParseDotBracketSecondary(pos.front(), pos.back());
  std::unique_ptr<energy::Structure> structure;
  printf("Energy: %d\n", energy::ComputeEnergy(secondary, em, &structure).energy);
  const auto descs = structure->Description();
  for (const auto& desc : descs) {
    printf("%s\n", desc.c_str());
  }
}
