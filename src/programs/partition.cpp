// Copyright 2016 E.
#include "partition/partition.h"

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>

#include "bridge/rnastructure.h"
#include "energy/load_model.h"
#include "fold/brute_fold.h"
#include "parsing.h"
#include "partition/partition_globals.h"

using namespace memerna;

void PrintProbabilities(const partition::probabilities_t& p) {
  const int N = static_cast<int>(p.Size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) { std::cout << std::setprecision(20) << p[i][j][0] << ' '; }
    std::cout << '\n';
  }
}

void PrintPartition(const partition::partition_t& p) {
  const int N = static_cast<int>(p.p.Size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) { std::cout << p.p[i][j][0] << ' '; }
    std::cout << '\n';
  }
}

int main(int argc, char* argv[]) {
  ArgParse argparse;
  argparse.AddOptions(energy::ENERGY_OPTIONS);
  argparse.AddOptions(CONTEXT_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  const auto& pos = argparse.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");
  const auto primary = parsing::StringToPrimary(pos.front());
  const auto em = energy::LoadEnergyModelFromArgParse(argparse);
  const bridge::Rnastructure rnastructure("extern/miles_rnastructure/data_tables/", false);

  Context ctx(primary, em, ContextOptionsFromArgParse(argparse));
  std::cout << "MEMERNA:\n";
  PrintProbabilities(partition::ComputeProbabilities(ctx.Partition()));
  std::cout << "RNAstructure:\n";
  PrintProbabilities(rnastructure.Partition(primary).second);
}
