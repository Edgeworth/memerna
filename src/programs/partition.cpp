// Copyright 2016 Eliot Courtney.
#include "compute/partition/partition.h"

#include <fmt/core.h>

#include "ctx/ctx.h"
#include "ctx/ctx_cfg.h"
#include "model/primary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  std::ios_base::sync_with_stdio(false);
  mrna::ArgParse args;
  mrna::ctx::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 1, "need primary sequence to fold");
  auto r = mrna::Primary::FromSeq(args.Pos(0));

  auto ctx = mrna::ctx::Ctx::FromArgParse(args);
  auto res = ctx.Partition(r);
  fmt::print("q: {}\np:\n", res.part.q);
  PrintPartition(res.part);
  fmt::print("\nprobabilities:\n");
  PrintBoltzProbs(res.prob);
}
