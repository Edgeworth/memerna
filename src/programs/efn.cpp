// Copyright 2016 Eliot Courtney.
#include <cstdio>

#include "compute/energy/load_model.h"
#include "compute/energy/structure.h"
#include "model/parsing.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args(mrna::energy::ENERGY_OPTS);
  args.AddOptions({{"v", {"verbose"}}});
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 2, "requires primary sequence and dot bracket");

  const auto em = mrna::energy::LoadEnergyModelFromArgParse(args);
  std::unique_ptr<mrna::energy::Structure> structure;
  mrna::computed_t computed;
  if (mrna::IsCtdString(pos.back())) {
    computed = mrna::ParseCtdComputed(pos.front(), pos.back());
    printf("Energy: %d\n",
        mrna::energy::ComputeEnergyWithCtds(computed, *em, false, &structure).energy);
  } else {
    const auto secondary = mrna::ParseDotBracketSecondary(pos.front(), pos.back());
    computed = mrna::energy::ComputeEnergy(secondary, *em, &structure);
    printf("Energy: %d\n", computed.energy);
  }

  if (args.HasFlag("v")) {
    printf("%s\n", mrna::ComputedToCtdString(computed).c_str());
    const auto descs = structure->Description();
    for (const auto& desc : descs) { printf("%s\n", desc.c_str()); }
  }
}
