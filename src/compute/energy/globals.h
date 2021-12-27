// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_GLOBALS_H_
#define COMPUTE_ENERGY_GLOBALS_H_

#include "compute/energy/fast_energy.h"
#include "compute/energy/model.h"

namespace mrna::energy {

extern EnergyModel gem;
extern precomp_t gpc;
void SetEnergyGlobalState(const primary_t& r, const EnergyModel& em);

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_GLOBALS_H_
