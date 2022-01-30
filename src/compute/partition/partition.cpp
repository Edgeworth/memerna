// Copyright 2016 Eliot Courtney.
#include "compute/partition/partition.h"

#include "compute/boltz_dp.h"

namespace mrna::part {

BoltzProbs Part::Prob() const {
  const int N = static_cast<int>(p.size());
  BoltzProbs prob(N, 0);
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j) prob[i][j] = p[i][j] * p[j][i] / q;
  return prob;
}

}  // namespace mrna::part
