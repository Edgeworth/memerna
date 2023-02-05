// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <ios>

#include "api/part.h"
#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "model/primary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  std::ios_base::sync_with_stdio(false);
  mrna::ArgParse args;
  mrna::api::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 1, "need primary sequence to fold");
  auto r = mrna::Primary::FromSeq(args.Pos(0));

  auto ctx = mrna::api::Ctx::FromArgParse(args);
  auto res = ctx.Partition(r);
  fmt::print("q: {}\np:\n", res.part.q);
  PrintPartition(res.part);
  fmt::print("\nprobabilities:\n");
  PrintBoltzProbs(res.prob);
}
