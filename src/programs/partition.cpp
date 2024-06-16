// Copyright 2016 E.
#include <fmt/core.h>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/pfn.h"
#include "model/pfn.h"
#include "model/primary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 1, "need primary sequence to fold");
  auto r = mrna::Primary::FromSeq(args.Pos(0));

  auto ctx = mrna::Ctx::FromArgParse(args);
  auto res = ctx.Pfn(r);
  fmt::print("q: {}\np:\n", res.pfn.q);
  PrintPfn(res.pfn.p);
  fmt::print("\nprobabilities:\n");
  PrintBoltzProbs(res.pfn.prob);
}
