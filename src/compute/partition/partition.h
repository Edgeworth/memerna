// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_PARTITION_PARTITION_H_
#define COMPUTE_PARTITION_PARTITION_H_

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "compute/energy/boltzmann_model.h"
#include "compute/energy/model.h"
#include "compute/energy/precomp.h"
#include "model/primary.h"
#include "util/array.h"

namespace mrna::partition {

enum : int8_t { PT_P, PT_U, PT_U2, PT_U_WC, PT_U_GU, PT_U_RCOAX, PT_SIZE };

enum : int8_t {
  PTEXT_R,
  PTEXT_L,
  PTEXT_R_WC,  // Must start with a branch not involved in an interaction that is Watson-Crick
  PTEXT_R_GU,  // Must start with a branch not involved in an interaction that is GU
  PTEXT_R_RCOAX,  // Must start with a branch, that branch is involved backwards in a RCOAX stack.
  PTEXT_L_WC,
  PTEXT_L_GU,
  PTEXT_L_LCOAX,
  PTEXT_SIZE
};

typedef Array3D<BoltzEnergy, 1> Probabilities;

struct Partition {
  Array3D<BoltzEnergy, 1> p;
  BoltzEnergy q;
};

Probabilities ComputeProbabilities(const Partition& partition);

namespace internal {

void Partition1(const Primary& r, const energy::BoltzEnergyModel& bem);
void Partition0(const Primary& r, const energy::EnergyModel& em);
void Exterior(const Primary& r, const energy::EnergyModel& em);

// Only works with [0, 2N).
inline int FastMod(int a, int m) {
  if (a >= m) return a - m;
  return a;
}

}  // namespace internal

}  // namespace mrna::partition

#endif  // COMPUTE_PARTITION_PARTITION_H_
