// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
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

void PrintProbabilities(const partition::partition_t& p) {
  auto probs = partition::ComputeProbabilities(p);
  const int N = int(p.p.Size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) { std::cout << std::setprecision(20) << probs[i][j][0] << ' '; }
    std::cout << '\n';
  }
}

void PrintPartition(const partition::partition_t& p) {
  const int N = int(p.p.Size());
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
  verify_expr(pos.size() == 1, "need primary sequence to fold");
  const auto primary = parsing::StringToPrimary(pos.front());
  const auto em = energy::LoadEnergyModelFromArgParse(argparse);
  const bridge::Rnastructure rnastructure("extern/miles_rnastructure/data_tables/", false);

  Context ctx(parsing::StringToPrimary(pos.front()), em, ContextOptionsFromArgParse(argparse));
  PrintProbabilities(ctx.Partition());
}
