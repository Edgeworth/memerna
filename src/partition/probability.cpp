// Copyright 2016 E.
#include <energy/energy_model.h>

#include "partition/partition.h"

namespace mrna {
namespace partition {

probabilities_t ComputeProbabilities(const partition_t& partition) {
  const int N = static_cast<int>(partition.p.Size());
  probabilities_t probabilities(std::size_t(N), 0);
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j)
      probabilities[i][j][0] = partition.p[i][j][0] * partition.p[j][i][0] / partition.q;
  return probabilities;
}

}  // namespace partition
}  // namespace mrna
