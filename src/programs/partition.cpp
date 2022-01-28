// Copyright 2016 Eliot Courtney.
#include "compute/partition/partition.h"

#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "compute/energy/energy.h"
#include "compute/energy/model.h"
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

  const auto& pos = args.positional();
  verify(pos.size() == 1, "need primary sequence to fold");
  auto r = mrna::Primary::FromString(pos.front());

  auto ctx = mrna::Ctx::FromArgParse(args);
  PrintBoltzProbs(ctx.Partition(std::move(r)).prob);
}
