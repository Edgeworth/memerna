// Copyright 2016 Eliot Courtney.
#include "api/subopt/subopt.h"

#include <fmt/core.h>

#include <string>

#include "api/ctx/ctx.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/options.h"
#include "api/subopt/subopt_cfg.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"
#include "util/error.h"

inline const auto OPT_CTD_OUTPUT =
    mrna::Opt(mrna::Opt::FLAG).LongName("ctd-output").Help("if we should output CTD data");

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_QUIET);
  args.RegisterOpt(OPT_CTD_OUTPUT);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 1, "need primary sequence to fold");

  auto ctx = mrna::Ctx::FromArgParse(args);
  auto mfe_alg = args.MaybeGet<mrna::MfeAlg>(mrna::OPT_MFE_ALG);
  auto subopt_alg = args.MaybeGet<mrna::SuboptAlg>(mrna::OPT_SUBOPT_ALG);
  const bool should_print = !args.GetOr(mrna::OPT_QUIET);
  const bool ctd_data = args.GetOr(OPT_CTD_OUTPUT);
  const auto subopt_cfg = mrna::subopt::SuboptCfg::FromArgParse(args);
  auto r = mrna::Primary::FromSeq(args.Pos(0));
  auto energy_cfg = mrna::erg::EnergyCfg::FromArgParse(args);
  auto pf = mrna::erg::PseudofreeCfg::FromArgParse(args);
  pf.Verify(r);

  mrna::subopt::SuboptCallback fn = [](const mrna::subopt::SuboptResult&) {};
  if (should_print) {
    if (ctd_data) {
      fn = [energy_cfg](const mrna::subopt::SuboptResult& c) {
        fmt::print("{} {}\n", c.energy, energy_cfg.ToCtdString(c.tb.s, c.tb.ctd));
      };
    } else {
      fn = [](const mrna::subopt::SuboptResult& c) {
        fmt::print("{} {}\n", c.energy, c.tb.s.ToDb());
      };
    }
  }
  int strucs = ctx.Subopt(r, mfe_alg, subopt_alg, energy_cfg, pf, fn, subopt_cfg);
  fmt::print("{} suboptimal structures\n", strucs);
}
