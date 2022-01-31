// Copyright 2016 E.
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "compute/energy/config.h"
#include "compute/energy/energy.h"
#include "compute/energy/model.h"
#include "compute/energy/structure.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "options.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::energy::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_VERBOSE);
  args.ParseOrExit(argc, argv);
  verify(args.PosSize() == 2, "requires primary sequence and dot bracket");

  const auto em = mrna::energy::EnergyModel::FromArgParse(args);
  const auto& rs = args.Pos(0);
  const auto& ss = args.Pos(1);
  mrna::energy::EnergyResult res;
  if (mrna::Ctds::IsCtdString(ss)) {
    const auto [r, s, ctd] = mrna::ParsePrimaryCtdString(rs, ss);
    res = em->TotalEnergy(r, s, &ctd, true);
    printf("Energy: %d\n", res.energy);
    printf("%s\n", res.ctd.ToString(s).c_str());
  } else {
    const auto [r, s] = mrna::ParsePrimaryDotBracket(rs, ss);
    res = em->TotalEnergy(r, s, nullptr, true);
    printf("Energy: %d\n", res.energy);
    printf("%s\n", res.ctd.ToString(s).c_str());
  }

  if (args.Has(mrna::OPT_VERBOSE)) {
    const auto descs = res.struc->Description();
    for (const auto& desc : descs) printf("%s\n", desc.c_str());
  }
}
