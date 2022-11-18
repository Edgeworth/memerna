// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_PARTITION_PARTITION_H_
#define COMPUTE_PARTITION_PARTITION_H_

#include <cassert>
#include <tuple>

#include "compute/boltz_dp.h"
#include "compute/energy/boltzmann_model.h"
#include "compute/energy/model.h"
#include "model/constants.h"
#include "model/primary.h"

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

std::tuple<BoltzDpArray, BoltzExtArray> Partition1(
    const Primary& r, const energy::BoltzEnergyModelPtr& bem);
std::tuple<BoltzDpArray, BoltzExtArray> Partition0(
    const Primary& r, const energy::EnergyModelPtr& em);
BoltzExtArray Exterior(const Primary& r, const energy::EnergyModel& em, const BoltzDpArray& dp);

// Only works with [0, 2N).
inline int FastMod(int a, int m) {
  assert(a < 2 * m);
  if (a >= m) return a - m;
  return a;
}

}  // namespace mrna::part

#endif  // COMPUTE_PARTITION_PARTITION_H_
