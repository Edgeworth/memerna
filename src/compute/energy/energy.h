// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_ENERGY_H_
#define COMPUTE_ENERGY_ENERGY_H_

#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <utility>

#include "compute/energy/model.h"
#include "model/structure.h"

namespace mrna::energy {

class Structure;

// If (st, en) is not paired, treated as an exterior loop.
energy_t ComputeSubstructureEnergy(computed_t& computed, bool compute_ctds, int st, int en,
    const EnergyModel& em, std::unique_ptr<Structure>* s = nullptr);
computed_t ComputeEnergy(
    const secondary_t& secondary, const EnergyModel& em, std::unique_ptr<Structure>* s = nullptr);
computed_t ComputeEnergyWithCtds(const computed_t& computed, const EnergyModel& em,
    bool compute_ctds = false, std::unique_ptr<Structure>* s = nullptr);

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_ENERGY_H_
