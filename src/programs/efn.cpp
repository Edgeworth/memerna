// Copyright 2016 Eliot Courtney.
#include <memory>
#include <string>
#include <vector>

#include "compute/energy/energy_cfg.h"
#include "compute/energy/model.h"
#include "compute/energy/structure.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "options.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  std::ios_base::sync_with_stdio(false);
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
    const auto [r, s, ctd] = mrna::ParseSeqCtdString(rs, ss);
    res = em->TotalEnergy(r, s, &ctd, true);
    std::cout << res.energy << '\n';
    std::cout << res.ctd.ToString(s) << '\n';
  } else {
    const auto [r, s] = mrna::ParseSeqDb(rs, ss);
    res = em->TotalEnergy(r, s, nullptr, true);
    std::cout << res.energy << '\n';
    std::cout << res.ctd.ToString(s) << '\n';
  }

  if (args.GetOr(mrna::OPT_VERBOSE)) {
    const auto descs = res.struc->Description();
    for (const auto& desc : descs) std::cout << desc << '\n';
  }
}
