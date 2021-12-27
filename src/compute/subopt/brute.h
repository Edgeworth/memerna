// Copyright 2021 Eliot Courtney.
#ifndef COMPUTE_SUBOPT_BRUTE_H_
#define COMPUTE_SUBOPT_BRUTE_H_

#include <utility>
#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"

namespace mrna::subopt {

std::vector<computed_t> SuboptimalBruteForce(
    const primary_t& r, const energy::EnergyModel& em, int max_structures_);

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_BRUTE_H_
