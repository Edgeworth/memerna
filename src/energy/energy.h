#ifndef MEMERNA_ENERGY_H
#define MEMERNA_ENERGY_H

#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <utility>
#include "argparse.h"
#include "common.h"
#include "energy/energy_model.h"

namespace memerna {
namespace energy {

class Structure;

computed_t ComputeEnergy(
    const secondary_t& secondary, const EnergyModel& em, std::unique_ptr<Structure>* s = nullptr);
computed_t ComputeEnergyWithCtds(const computed_t& computed, const EnergyModel& em,
    bool compute_ctds = false, std::unique_ptr<Structure>* s = nullptr);
}
}

#endif  // MEMERNA_ENERGY_H
