// Copyright 2016 Eliot Courtney.
#include <cstdio>

#include "compute/energy/model.h"
#include "compute/energy/structure.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args(mrna::energy::ENERGY_OPTS);
  args.AddOptions({{"v", {"verbose"}}});
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 2, "requires primary sequence and dot bracket");

  const auto em = mrna::energy::EnergyModel::FromArgParse(args);
  mrna::energy::EnergyResult res;
  mrna::Secondary s;
  if (mrna::Ctds::IsCtdString(pos.back())) {
    const auto [r, s, ctd] = mrna::ParsePrimaryCtdString(pos.front(), pos.back());
    res = em.TotalEnergy(r, s, &ctd, true);
  } else {
    const auto [r, s] = mrna::ParsePrimaryDotBracket(pos.front(), pos.back());
    res = em.TotalEnergy(r, s, nullptr, true);
  }
  printf("Energy: %d\n", res.energy);

  if (args.HasFlag("v")) {
    printf("%s\n", res.ctd.ToString(s).c_str());
    const auto descs = res.struc->Description();
    for (const auto& desc : descs) printf("%s\n", desc.c_str());
  }
}
