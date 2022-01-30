// Copyright 2016 E.
#ifndef COMPUTE_PARTITION_PARTITION_H_
#define COMPUTE_PARTITION_PARTITION_H_

#include <cassert>
#include <tuple>

#include "compute/boltz_dp.h"
#include "compute/energy/boltzmann_model.h"
#include "compute/energy/model.h"
#include "model/model.h"
#include "model/primary.h"

namespace mrna::partition {

struct Partition {
  BoltzSums p;
  BoltzEnergy q;
};

struct PartitionResult {
  BoltzDpArray dp;
  BoltzExtArray ext;
  Partition p;
  BoltzProbs prob;
};

BoltzProbs ComputeBoltzProbs(const Partition& partition);

std::tuple<BoltzDpArray, BoltzExtArray> Partition1(
    const Primary& r, const energy::BoltzEnergyModel& bem);
std::tuple<BoltzDpArray, BoltzExtArray> Partition0(const Primary& r, const energy::EnergyModel& em);
BoltzExtArray Exterior(const Primary& r, const energy::EnergyModel& em, const BoltzDpArray& dp);

// Only works with [0, 2N).
inline int FastMod(int a, int m) {
  assert(a < 2 * m);
  if (a >= m) return a - m;
  return a;
}

}  // namespace mrna::partition

#endif  // COMPUTE_PARTITION_PARTITION_H_
