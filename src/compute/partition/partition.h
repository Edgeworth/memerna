// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_PARTITION_PARTITION_H_
#define COMPUTE_PARTITION_PARTITION_H_

#include <cassert>

#include "compute/boltz_dp.h"
#include "model/energy.h"

namespace mrna::part {

struct Part {
  BoltzSums p;
  BoltzEnergy q;

  [[nodiscard]] BoltzProbs Prob() const;
};

struct PartResult {
  BoltzDpArray dp;
  BoltzExtArray ext;
  Part part;
  BoltzProbs prob;
};

// Only works with [0, 2N).
inline int FastMod(int a, int m) {
  assert(a < 2 * m);
  if (a >= m) return a - m;
  return a;
}

}  // namespace mrna::part

#endif  // COMPUTE_PARTITION_PARTITION_H_
