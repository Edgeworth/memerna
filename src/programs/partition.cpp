// Copyright 2016 E.
#include "compute/partition/partition.h"

#include <iostream>
#include <utility>

#include "ctx/config.h"
#include "ctx/ctx.h"
#include "model/primary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::ctx::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 1, "need primary sequence to fold");
  auto r = mrna::Primary::FromString(args.Pos(0));

  auto ctx = mrna::ctx::Ctx::FromArgParse(args);
  auto res = ctx.Partition(std::move(r));
  std::cout << "q: " << res.part.q << '\n';
  std::cout << "p:\n";
  PrintPartition(res.part);
  std::cout << "\nprobabilities:\n";
  PrintBoltzProbs(res.prob);
}
