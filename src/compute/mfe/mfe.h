// Copyright 2016 E.
#ifndef COMPUTE_MFE_MFE_H_
#define COMPUTE_MFE_MFE_H_

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "model/model.h"
#include "model/primary.h"

namespace mrna::mfe {

struct MfeResult {
  DpArray dp;
  ExtArray ext;
  Energy energy = 0;
};

// TODO(3): Move this?
struct Cand {
  Energy energy;
  int idx;
};

// Basic MFE folding.
DpArray ComputeTables0(const Primary& r, const energy::EnergyModelPtr& em);

// Basic MFE folding.
DpArray ComputeTables1(const Primary& r, const energy::EnergyModelPtr& em);

// Sparse folding.
DpArray ComputeTables2(const Primary& r, const energy::EnergyModelPtr& em);

// Sparse folding with Lyngso's algorithm.
DpArray ComputeTables3(const Primary& r, const energy::EnergyModelPtr& em);
ExtArray ComputeExterior(const Primary& r, const energy::EnergyModel& em, const DpArray& dp);

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_MFE_H_
