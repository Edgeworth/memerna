// Copyright 2016 E.
#include <cstdio>

#include "compute/energy/energy.h"
#include "compute/energy/load_model.h"
#include "compute/energy/structure.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args(mrna::energy::ENERGY_OPTS);
  args.AddOptions({{"v", {"verbose"}}});
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 2, "requires primary sequence and dot bracket");

  const auto em = mrna::energy::LoadEnergyModelFromArgParse(args);
  std::unique_ptr<mrna::energy::Structure> struc;
  mrna::energy::EnergyResult res;
  mrna::Secondary s;
  if (mrna::IsCtdString(pos.back())) {
    const auto [r, s, ctd] = mrna::ParsePrimaryCtdString(pos.front(), pos.back());
    res = mrna::energy::ComputeEnergy(r, s, &ctd, em, &struc);
  } else {
    const auto [r, s] = mrna::ParsePrimaryDotBracket(pos.front(), pos.back());
    res = mrna::energy::ComputeEnergy(r, s, nullptr, em, &struc);
  }
  printf("Energy: %d\n", res.energy);

  if (args.HasFlag("v")) {
    printf("%s\n", mrna::CtdString(s, res.ctd).c_str());
    const auto descs = struc->Description();
    for (const auto& desc : descs) { printf("%s\n", desc.c_str()); }
  }
}
