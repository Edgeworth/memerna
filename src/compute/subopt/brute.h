// Copyright 2021 E.
#ifndef COMPUTE_SUBOPT_BRUTE_H_
#define COMPUTE_SUBOPT_BRUTE_H_

#include <utility>
#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"

namespace mrna::subopt {

std::vector<subopt::SuboptResult> SuboptimalBruteForce(
    Primary r, const energy::EnergyModel& em, int max_structures);

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_BRUTE_H_
