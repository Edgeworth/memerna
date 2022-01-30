// Copyright 2016 E.
#include "compute/subopt/subopt.h"

#include <cstdio>
#include <string>

#include "ctx/config.h"
#include "ctx/ctx.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "options.h"
#include "util/argparse.h"
#include "util/error.h"
#include "compute/subopt/config.h"

// TODO: can replace this with modelcfg ctd stuff
inline const auto OPT_CTD_OUTPUT =
    mrna::Opt().LongName("ctd-output").Help("if we should output CTD data");

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::ctx::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_QUIET);
  args.RegisterOpt(OPT_CTD_OUTPUT);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 1, "need primary sequence to fold");

  auto ctx = mrna::ctx::Ctx::FromArgParse(args);
  const bool should_print = !args.Has(mrna::OPT_QUIET);
  const bool ctd_data = args.Has(OPT_CTD_OUTPUT);
  const auto cfg = mrna::subopt::SuboptCfg::FromArgParse(args);

  mrna::subopt::SuboptCallback fn = [](const mrna::subopt::SuboptResult&) {};
  if (should_print) {
    if (ctd_data) {
      fn = [](const mrna::subopt::SuboptResult& c) {
        printf("%d ", c.energy);
        puts(c.tb.ctd.ToString(c.tb.s).c_str());
      };
    } else {
      fn = [](const mrna::subopt::SuboptResult& c) {
        printf("%d ", c.energy);
        puts(c.tb.s.ToDotBracket().c_str());
      };
    }
  }
  int strucs = ctx.Suboptimal(mrna::Primary::FromString(args.Pos(0)), fn, cfg);
  printf("%d suboptimal structures\n", strucs);
}
