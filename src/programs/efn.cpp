// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <string>

#include "api/ctx/ctx.h"
#include "api/energy/pseudofree_cfg.h"
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
  mrna::RegisterOpts(&args);
  args.RegisterOpt(OPT_DETAIL);
  args.ParseOrExit(argc, argv);
  verify(args.PosSize() == 2, "requires primary sequence and dot bracket");

  auto ctx = mrna::Ctx::FromArgParse(args);
  const auto energy_cfg = mrna::erg::EnergyCfg::FromArgParse(args);
  auto pf = mrna::erg::PseudofreeCfg::FromArgParse(args);
  const auto& rs = args.Pos(0);
  const auto& ss = args.Pos(1);
  mrna::erg::EnergyResult res;
  if (mrna::Ctds::IsCtdString(ss)) {
    const auto [r, s, ctd] = energy_cfg.ParseSeqCtdString(rs, ss);
    pf.Verify(r);
    res = ctx.Efn(r, s, energy_cfg, pf, &ctd, /*build_structure=*/true);
    fmt::print("{}\n", res.energy);
    fmt::print("{}\n", energy_cfg.ToCtdString(s, res.ctd));
  } else {
    const auto [r, s] = mrna::ParseSeqDb(rs, ss);
    pf.Verify(r);
    res = ctx.Efn(r, s, energy_cfg, pf, /*given_ctd=*/nullptr, /*build_structure=*/true);
    fmt::print("{}\n", res.energy);
    fmt::print("{}\n", energy_cfg.ToCtdString(s, res.ctd));
  }

  if (args.GetOr(OPT_DETAIL)) {
    const auto descs = res.struc->Description();
    for (const auto& desc : descs) fmt::print("{}\n", desc);
  }
}
