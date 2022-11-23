// Copyright 2016 Eliot Courtney.
#include "compute/subopt/subopt.h"

#include <iostream>
#include <string>

#include "compute/subopt/subopt_cfg.h"
#include "ctx/ctx.h"
#include "ctx/ctx_cfg.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "options.h"
#include "util/argparse.h"
#include "util/error.h"

// TODO(2): can replace this with modelcfg ctd stuff
inline const auto OPT_CTD_OUTPUT =
    mrna::Opt(mrna::Opt::FLAG).LongName("ctd-output").Help("if we should output CTD data");

int main(int argc, char* argv[]) {
  std::ios_base::sync_with_stdio(false);
  mrna::ArgParse args;
  mrna::ctx::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_QUIET);
  args.RegisterOpt(OPT_CTD_OUTPUT);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 1, "need primary sequence to fold");

  auto ctx = mrna::ctx::Ctx::FromArgParse(args);
  const bool should_print = !args.GetOr(mrna::OPT_QUIET);
  const bool ctd_data = args.GetOr(OPT_CTD_OUTPUT);
  const auto cfg = mrna::subopt::SuboptCfg::FromArgParse(args);

  mrna::subopt::SuboptCallback fn = [](const mrna::subopt::SuboptResult&) {};
  if (should_print) {
    if (ctd_data) {
      fn = [](const mrna::subopt::SuboptResult& c) {
        std::cout << c.energy << ' ' << c.tb.ctd.ToString(c.tb.s) << '\n';
      };
    } else {
      fn = [](const mrna::subopt::SuboptResult& c) {
        std::cout << c.energy << ' ' << c.tb.s.ToDb() << '\n';
      };
    }
  }
  int strucs = ctx.Suboptimal(mrna::Primary::FromSeq(args.Pos(0)), fn, cfg);
  std::cout << strucs << " suboptimal structures\n";
}
