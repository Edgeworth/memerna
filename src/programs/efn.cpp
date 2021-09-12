// Copyright 2016 E.
#include <cstdio>

#include "energy/load_model.h"
#include "energy/structure.h"
#include "parsing.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions({{"v", {"verbose"}}});
  argparse.ParseOrExit(argc, argv);
  const auto& pos = argparse.GetPositional();
  verify_expr(pos.size() == 2, "requires primary sequence and dot bracket");

  const auto em = energy::LoadEnergyModelFromArgParse(argparse);
  std::unique_ptr<energy::Structure> structure;
  computed_t computed;
  if (parsing::IsCtdString(pos.back())) {
    computed = parsing::ParseCtdComputed(pos.front(), pos.back());
    printf("Energy: %d\n", energy::ComputeEnergyWithCtds(computed, *em, false, &structure).energy);
  } else {
    const auto secondary = parsing::ParseDotBracketSecondary(pos.front(), pos.back());
    computed = energy::ComputeEnergy(secondary, *em, &structure);
    printf("Energy: %d\n", computed.energy);
  }

  if (argparse.HasFlag("v")) {
    printf("%s\n", parsing::ComputedToCtdString(computed).c_str());
    const auto descs = structure->Description();
    for (const auto& desc : descs) { printf("%s\n", desc.c_str()); }
  }
}
