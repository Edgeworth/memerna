// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <memory>
#include <string>

#include "api/ctx/backend.h"
#include "api/ctx/backend_cfg.h"
#include "api/energy/energy.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "model/structure.h"
#include "util/argparse.h"
#include "util/error.h"

inline const mrna::Opt OPT_DETAIL =
    mrna::Opt(mrna::Opt::FLAG).LongName("detail").ShortName("d").Help("detailed structure output");

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::RegisterOptsBackendCfg(&args);
  args.RegisterOpt(OPT_DETAIL);
  args.ParseOrExit(argc, argv);
  verify(args.PosSize() == 2, "requires primary sequence and dot bracket");

  const auto m = mrna::BackendFromArgParse(args);
  const auto& rs = args.Pos(0);
  const auto& ss = args.Pos(1);
  mrna::erg::EnergyResult res;
  if (mrna::Ctds::IsCtdString(ss)) {
    const auto [r, s, ctd] = mrna::ParseSeqCtdString(rs, ss);
    res = mrna::TotalEnergy(m, r, s, &ctd, true);
    fmt::print("{}\n", res.energy);
    fmt::print("{}\n", res.ctd.ToString(s));
  } else {
    const auto [r, s] = mrna::ParseSeqDb(rs, ss);
    res = mrna::TotalEnergy(m, r, s, nullptr, true);
    fmt::print("{}\n", res.energy);
    fmt::print("{}\n", res.ctd.ToString(s));
  }

  if (args.GetOr(OPT_DETAIL)) {
    const auto descs = res.struc->Description();
    for (const auto& desc : descs) fmt::print("{}\n", desc);
  }
}
