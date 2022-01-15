// Copyright 2016 Eliot Courtney.
#include "compute/partition/partition.h"

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>

#include "compute/energy/load_model.h"
#include "compute/mfe/brute.h"
#include "model/context.h"

void PrintProbabilities(const mrna::Probabilities& p) {
  const int N = static_cast<int>(p.Size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) { std::cout << std::setprecision(20) << p[i][j][0] << ' '; }
    std::cout << '\n';
  }
}

void PrintPartition(const mrna::partition::Partition& p) {
  const int N = static_cast<int>(p.p.Size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) { std::cout << p.p[i][j][0] << ' '; }
    std::cout << '\n';
  }
}

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  args.AddOptions(mrna::energy::ENERGY_OPTS);
  args.AddOptions(mrna::MODEL_OPTS);
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");
  auto r = mrna::Primary::FromString(pos.front());
  const auto em = mrna::energy::LoadEnergyModelFromArgParse(args);

  mrna::Context ctx(std::move(r), em, ModelCfgFromArgParse(args));
  PrintProbabilities(ctx.Partition().prob);
}
