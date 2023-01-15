// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <memory>
#include <string>
#include <vector>

#include "compute/energy/energy.h"
#include "compute/energy/energy_cfg.h"
#include "compute/energy/model.h"
#include "compute/energy/structure.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/secondary.h"
#include "options.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  std::ios_base::sync_with_stdio(false);
  mrna::ArgParse args;
  mrna::erg::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_VERBOSE);
  args.ParseOrExit(argc, argv);
  verify(args.PosSize() == 2, "requires primary sequence and dot bracket");

  const auto em = mrna::erg::FromArgParse(args);
  const auto& rs = args.Pos(0);
  const auto& ss = args.Pos(1);
  mrna::erg::EnergyResult res;
  if (mrna::Ctds::IsCtdString(ss)) {
    const auto [r, s, ctd] = mrna::ParseSeqCtdString(rs, ss);
    res = mrna::erg::TotalEnergy(em, r, s, &ctd, true);
    fmt::print("{}\n", res.energy);
    fmt::print("{}\n", res.ctd.ToString(s));
  } else {
    const auto [r, s] = mrna::ParseSeqDb(rs, ss);
    res = mrna::erg::TotalEnergy(em, r, s, nullptr, true);
    fmt::print("{}\n", res.energy);
    fmt::print("{}\n", res.ctd.ToString(s));
  }

  if (args.GetOr(mrna::OPT_VERBOSE)) {
    const auto descs = res.struc->Description();
    for (const auto& desc : descs) fmt::print("{}\n", desc);
  }
}
