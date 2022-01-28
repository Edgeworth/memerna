// Copyright 2016 Eliot Courtney.
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "compute/energy/energy.h"
#include "compute/energy/model.h"
#include "compute/energy/structure.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::energy::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_VERBOSE);
  args.ParseOrExit(argc, argv);
  const auto& pos = args.positional();
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

  if (args.Has(mrna::OPT_VERBOSE)) {
    printf("%s\n", res.ctd.ToString(s).c_str());
    const auto descs = res.struc->Description();
    for (const auto& desc : descs) printf("%s\n", desc.c_str());
  }
}
