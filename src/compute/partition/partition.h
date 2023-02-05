// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_PARTITION_PARTITION_H_
#define COMPUTE_PARTITION_PARTITION_H_

#include <cassert>
#include <variant>

#include "compute/partition/t04/dp.h"
#include "model/energy.h"

namespace mrna::part {

// Holds the Boltzmann sums used in the partition function.
using BoltzSums = Array2D<BoltzEnergy>;

using BoltzProbs = Array2D<BoltzEnergy>;

using DpState = std::variant<t04::DpState>;

struct Part {
  BoltzSums p;
  BoltzEnergy q;

  [[nodiscard]] BoltzProbs Prob() const;
};

struct PartResult {
  DpState dp;
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
