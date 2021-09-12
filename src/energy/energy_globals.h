// Copyright 2016 E.
#ifndef MEMERNA_ENERGY_GLOBALS_H_
#define MEMERNA_ENERGY_GLOBALS_H_

#include "common.h"
#include "energy/energy_model.h"
#include "energy/fast_energy.h"

namespace memerna {
namespace energy {

extern EnergyModel gem;
extern precomp_t gpc;
void SetEnergyGlobalState(const primary_t& r, const EnergyModel& em);

}  // namespace energy
}  // namespace memerna

#endif  // MEMERNA_ENERGY_GLOBALS_H_
