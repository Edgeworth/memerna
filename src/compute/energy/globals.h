// Copyright 2016 E.
#ifndef COMPUTE_ENERGY_GLOBALS_H_
#define COMPUTE_ENERGY_GLOBALS_H_

#include "compute/energy/fast_energy.h"
#include "compute/energy/model.h"

namespace mrna::energy {

extern EnergyModel gem;
extern Precomp gpc;
void SetEnergyGlobalState(const Primary& r, const EnergyModel& em);

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_GLOBALS_H_
