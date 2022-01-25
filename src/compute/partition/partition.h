// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_PARTITION_PARTITION_H_
#define COMPUTE_PARTITION_PARTITION_H_

#include <cassert>
#include <tuple>

#include "compute/boltz_dp.h"
#include "compute/energy/boltzmann_model.h"
#include "model/model.h"
#include "model/primary.h"
#include "util/array.h"

namespace mrna::partition {

struct Partition {
  Array3D<BoltzEnergy, 1> p;
  BoltzEnergy q;
};

struct PartitionResult {
  Partition p;
  Probabilities prob;
};

Probabilities ComputeProbabilities(const Partition& partition);

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
