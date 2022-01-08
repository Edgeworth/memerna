// Copyright 2016 E.
#ifndef COMPUTE_MFE_MFE_H_
#define COMPUTE_MFE_MFE_H_

#include "compute/mfe/globals.h"
#include "model/base.h"

namespace mrna::mfe::internal {

struct Cand {
  Energy energy;
  int idx;
};

void ComputeTables0(const Primary& r, const energy::EnergyModel& em);
void ComputeTables1(const Primary& r, const energy::EnergyModel& em);
void ComputeTables2(const Primary& r, const energy::EnergyModel& em);
void ComputeTables3(const Primary& r, const energy::EnergyModel& em);
void ComputeExterior(const Primary& r, const energy::EnergyModel& em);

}  // namespace mrna::mfe::internal

#endif  // COMPUTE_MFE_MFE_H_
