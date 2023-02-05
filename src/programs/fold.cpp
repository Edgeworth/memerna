// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <algorithm>
#include <ios>
#include <string>

#include "api/mfe.h"
#include "api/trace.h"
#include "ctx/ctx.h"
#include "ctx/ctx_cfg.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  std::ios_base::sync_with_stdio(false);
  mrna::ArgParse args;
  mrna::api::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);
  verify(args.PosSize() == 1, "need primary sequence to fold");

  auto ctx = mrna::api::Ctx::FromArgParse(args);
  const auto res = ctx.Fold(mrna::Primary::FromSeq(args.Pos(0)));

  fmt::print("{}\n", res.mfe.energy);
  fmt::print("{}\n", res.tb.s.ToDb());
  fmt::print("{}\n", res.tb.ctd.ToString(res.tb.s));
}
