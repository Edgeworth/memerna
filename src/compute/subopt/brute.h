// Copyright 2021 E.
#ifndef COMPUTE_SUBOPT_BRUTE_H_
#define COMPUTE_SUBOPT_BRUTE_H_

#include <utility>
#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"

namespace mrna::subopt {

std::vector<Computed> SuboptimalBruteForce(
    const Primary& r, const energy::EnergyModel& em, int max_structures_);

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_BRUTE_H_
