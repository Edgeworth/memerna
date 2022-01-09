// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_ENERGY_H_
#define COMPUTE_ENERGY_ENERGY_H_

#include <cmath>
#include <deque>
#include <memory>
#include <optional>
#include <string>
#include <utility>

#include "compute/energy/model.h"
#include "model/ctd.h"
#include "model/model.h"

namespace mrna::energy {

struct EnergyResult {
  Energy energy = 0;
  // TODO: tree, remove struc args
  Ctds ctd;  // May be empty if CTDs were not computed.
};

class Structure;

// If (st, en) is not paired, treated as an exterior loop.
// If |ctd| is non-null, use the given ctds.
EnergyResult ComputeSubstructureEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
    int st, int en, const EnergyModel& em, std::unique_ptr<Structure>* struc = nullptr);
EnergyResult ComputeEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
    const EnergyModel& em, std::unique_ptr<Structure>* struc = nullptr);
}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_ENERGY_H_
