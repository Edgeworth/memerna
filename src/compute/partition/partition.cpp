// Copyright 2016 E.
#include "compute/partition/partition.h"

#include "compute/boltz_dp.h"

namespace mrna::partition {

BoltzProbs ComputeBoltzProbs(const Partition& p) {
  const int N = static_cast<int>(p.p.size());
  BoltzProbs probs(N, 0);
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j) probs[i][j] = p.p[i][j] * p.p[j][i] / p.q;
  return probs;
}

}  // namespace mrna::partition
