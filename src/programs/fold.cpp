// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/mfe.h"
#include "api/trace/trace.h"
#include "api/trace/trace_cfg.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);
  verify(args.PosSize() == 1, "need primary sequence to fold");

  auto ctx = mrna::Ctx::FromArgParse(args);
  auto trace_cfg = mrna::trace::TraceCfg::FromArgParse(args);
  const auto res = ctx.Fold(mrna::Primary::FromSeq(args.Pos(0)), trace_cfg);

  fmt::print("{}\n", res.mfe.energy);
  fmt::print("{}\n", res.tb.s.ToDb());
  fmt::print("{}\n", mrna::BackendEnergyCfg(ctx.m()).ToCtdString(res.tb.s, res.tb.ctd));
}
