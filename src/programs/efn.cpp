// Copyright 2016 E.
#include <cstdio>

#include "energy/load_model.h"
#include "energy/structure.h"
#include "parsing.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse argparse(mrna::energy::ENERGY_OPTIONS);
  argparse.AddOptions({{"v", {"verbose"}}});
  argparse.ParseOrExit(argc, argv);
  const auto& pos = argparse.GetPositional();
  verify(pos.size() == 2, "requires primary sequence and dot bracket");

  const auto em = mrna::energy::LoadEnergyModelFromArgParse(argparse);
  std::unique_ptr<mrna::energy::Structure> structure;
  mrna::computed_t computed;
  if (mrna::parsing::IsCtdString(pos.back())) {
    computed = mrna::parsing::ParseCtdComputed(pos.front(), pos.back());
    printf("Energy: %d\n",
        mrna::energy::ComputeEnergyWithCtds(computed, *em, false, &structure).energy);
  } else {
    const auto secondary = mrna::parsing::ParseDotBracketSecondary(pos.front(), pos.back());
    computed = mrna::energy::ComputeEnergy(secondary, *em, &structure);
    printf("Energy: %d\n", computed.energy);
  }

  if (argparse.HasFlag("v")) {
    printf("%s\n", mrna::parsing::ComputedToCtdString(computed).c_str());
    const auto descs = structure->Description();
    for (const auto& desc : descs) { printf("%s\n", desc.c_str()); }
  }
}
