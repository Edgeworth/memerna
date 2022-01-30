// Copyright 2022 Eliot Courtney.
#include "programs/print.h"

#include <iomanip>
#include <iostream>

void PrintBoltzProbs(const mrna::BoltzProbs& p) {
  const int N = static_cast<int>(p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) std::cout << std::setprecision(20) << p[i][j] << ' ';
    std::cout << '\n';
  }
}
void PrintPartition(const mrna::part::Part& p) {
  const int N = static_cast<int>(p.p.size());
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) std::cout << p.p[i][j] << ' ';
    std::cout << '\n';
  }
}
