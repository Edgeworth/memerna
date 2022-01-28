// Copyright 2016 Eliot Courtney.
#include "compute/subopt/subopt.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "compute/energy/energy.h"
#include "compute/energy/model.h"
#include "context/config.h"
#include "context/ctx.h"
#include "fuzz/config.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"
#include "util/error.h"

inline const auto OPT_DELTA =
    mrna::Opt().LongName("delta").Default("-1").Help("maximum energy delta from minimum");
inline const auto OPT_NUM =
    mrna::Opt().LongName("num").Default("-1").Help("maximum number of reported structures");
inline const auto OPT_SORTED =
    mrna::Opt().LongName("sorted").Help("if the structures should be sorted");
inline const auto OPT_CTD_OUTPUT =
    mrna::Opt().LongName("ctd-output").Help("if we should output CTD data");

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::fuzz::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_QUIET);
  args.RegisterOpt(OPT_DELTA);
  args.RegisterOpt(OPT_NUM);
  args.RegisterOpt(OPT_SORTED);
  args.RegisterOpt(OPT_CTD_OUTPUT);
  args.ParseOrExit(argc, argv);

  const auto& pos = args.positional();
  verify(pos.size() == 1, "need primary sequence to fold");

  auto ctx = mrna::Ctx::FromArgParse(args);

  const mrna::Energy subopt_delta = atoi(args.Get(OPT_DELTA).c_str());
  const int subopt_num = atoi(args.Get(OPT_NUM).c_str());
  const bool should_print = !args.Has(mrna::OPT_QUIET);
  const bool sorted = args.Has(OPT_SORTED);
  const bool ctd_data = args.Has(OPT_CTD_OUTPUT);
  verify(subopt_delta >= 0 || subopt_num > 0, "nothing to do");

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
  int num_structures =
      ctx.Suboptimal(mrna::Primary::FromString(pos.front()), fn, sorted, subopt_delta, subopt_num);
  printf("%d suboptimal structures\n", num_structures);
}
