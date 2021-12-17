// Copyright 2016 E.
#include "compute/partition/partition.h"

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>

#include "bridge/rnastructure.h"
#include "compute/energy/load_model.h"
#include "compute/mfe/brute.h"
#include "compute/partition/globals.h"
#include "model/parsing.h"

void PrintProbabilities(const mrna::partition::probabilities_t& p) {
  const int N = static_cast<int>(p.Size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) { std::cout << std::setprecision(20) << p[i][j][0] << ' '; }
    std::cout << '\n';
  }
}

void PrintPartition(const mrna::partition::partition_t& p) {
  const int N = static_cast<int>(p.p.Size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) { std::cout << p.p[i][j][0] << ' '; }
    std::cout << '\n';
  }
}

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  args.AddOptions(mrna::energy::COMPUTE_ENERGY_OPTIONS);
  args.AddOptions(mrna::CONTEXT_OPTIONS);
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");
  const auto primary = mrna::StringToPrimary(pos.front());
  const auto em = mrna::energy::LoadEnergyModelFromArgParse(args);
  const mrna::bridge::RNAstructure rnastructure(args.GetOption("rnastructure-data"), false);

  mrna::Context ctx(primary, em, ContextOptionsFromArgParse(args));
  std::cout << "MEMERNA:\n";
  PrintProbabilities(mrna::partition::ComputeProbabilities(ctx.Partition()));
  std::cout << "RNAstructure:\n";
  PrintProbabilities(rnastructure.Partition(primary).second);
}
