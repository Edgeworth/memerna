// Copyright 2016 Eliot Courtney.
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

// TODO: Move this?
struct Cand {
  Energy energy;
  int idx;
};

DpArray ComputeTables0(const Primary& r, const energy::EnergyModel& em);
DpArray ComputeTables1(const Primary& r, const energy::EnergyModel& em);
DpArray ComputeTables2(const Primary& r, const energy::EnergyModel& em);
DpArray ComputeTables3(const Primary& r, const energy::EnergyModel& em);
ExtArray ComputeExterior(const Primary& r, const energy::EnergyModel& em, const DpArray& dp);

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_MFE_H_
