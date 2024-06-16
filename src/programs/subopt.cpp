// Copyright 2016 Eliot Courtney.
#include "api/subopt/subopt.h"

#include <fmt/core.h>

#include <string>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/options.h"
#include "api/subopt/subopt_cfg.h"
#include "model/ctd.h"
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
  const bool should_print = !args.GetOr(mrna::OPT_QUIET);
  const bool ctd_data = args.GetOr(OPT_CTD_OUTPUT);
  const auto cfg = mrna::subopt::SuboptCfg::FromArgParse(args);

  mrna::subopt::SuboptCallback fn = [](const mrna::subopt::SuboptResult&) {};
  if (should_print) {
    if (ctd_data) {
      fn = [](const mrna::subopt::SuboptResult& c) {
        fmt::print("{} {}\n", c.energy, c.tb.ctd.ToString(c.tb.s));
      };
    } else {
      fn = [](const mrna::subopt::SuboptResult& c) {
        fmt::print("{} {}\n", c.energy, c.tb.s.ToDb());
      };
    }
  }
  int strucs = ctx.Subopt(mrna::Primary::FromSeq(args.Pos(0)), fn, cfg);
  fmt::print("{} suboptimal structures\n", strucs);
}
