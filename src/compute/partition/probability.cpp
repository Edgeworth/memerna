// Copyright 2016 E.
#include <cstddef>

#include "compute/boltz_dp.h"
#include "compute/partition/partition.h"
#include "model/model.h"
#include "util/array.h"

namespace mrna::partition {

Probabilities ComputeProbabilities(const Partition& partition) {
  const int N = static_cast<int>(partition.p.size());
  Probabilities probabilities(std::size_t(N), 0);
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j)
      probabilities[i][j][0] = partition.p[i][j][0] * partition.p[j][i][0] / partition.q;
  return probabilities;
}

}  // namespace mrna::partition
