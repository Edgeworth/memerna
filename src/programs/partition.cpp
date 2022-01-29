// Copyright 2016 E.
#include "compute/partition/partition.h"

#include <iomanip>
#include <iostream>
#include <utility>

#include "context/config.h"
#include "context/ctx.h"
#include "model/primary.h"
#include "util/argparse.h"
#include "util/error.h"

void PrintBoltzProbs(const mrna::BoltzProbs& p) {
  const int N = static_cast<int>(p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) std::cout << std::setprecision(20) << p[i][j] << ' ';
    std::cout << '\n';
  }
}

void PrintPartition(const mrna::partition::Partition& p) {
  const int N = static_cast<int>(p.p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) std::cout << p.p[i][j] << ' ';
    std::cout << '\n';
  }
}

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 1, "need primary sequence to fold");
  auto r = mrna::Primary::FromString(args.Pos(0));

  auto ctx = mrna::Ctx::FromArgParse(args);
  auto res = ctx.Partition(std::move(r));
  std::cout << "q: " << res.p.q << '\n';
  std::cout << "p:\n";
  PrintPartition(res.p);
  std::cout << "\nprobabilities:\n";
  PrintBoltzProbs(res.prob);
}
