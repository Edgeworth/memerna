// Copyright 2016 E.
#ifndef COMPUTE_ENERGY_ENERGY_H_
#define COMPUTE_ENERGY_ENERGY_H_

#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <utility>

#include "compute/energy/model.h"
#include "model/ctd.h"
#include "model/model.h"

namespace mrna::energy {

class Structure;

// If (st, en) is not paired, treated as an exterior loop.
Energy ComputeSubstructureEnergy(Computed& computed, bool compute_ctds, int st, int en,
    const EnergyModel& em, std::unique_ptr<Structure>* s = nullptr);
Computed ComputeEnergy(
    const Secondary& secondary, const EnergyModel& em, std::unique_ptr<Structure>* s = nullptr);
Computed ComputeEnergyWithCtds(const Computed& computed, const EnergyModel& em,
    bool compute_ctds = false, std::unique_ptr<Structure>* s = nullptr);

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_ENERGY_H_
