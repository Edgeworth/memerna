// Copyright 2022 Eliot Courtney.
#include "programs/print.h"

#include <fmt/core.h>

void PrintBoltzProbs(const mrna::BoltzProbs& p) {
  const int N = static_cast<int>(p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      fmt::print("{:.20f} ", p[i][j]);
    }
    fmt::print("\n");
  }
}
void PrintPartition(const mrna::part::Part& p) {
  const int N = static_cast<int>(p.p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      fmt::print("{} ", p.p[i][j]);
    }
    fmt::print("\n");
  }
}
