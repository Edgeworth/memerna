// Copyright 2022 E.
#ifndef COMPUTE_BRUTE_ALG_H_
#define COMPUTE_BRUTE_ALG_H_

#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/config.h"
#include "compute/subopt/subopt.h"
#include "model/primary.h"

namespace mrna::brute {

subopt::SuboptResult MfeBruteForce(Primary r, energy::EnergyModelPtr em);

part::PartResult PartitionBruteForce(Primary r, energy::EnergyModelPtr em);

std::vector<subopt::SuboptResult> SuboptimalBruteForce(
    Primary r, energy::EnergyModelPtr em, subopt::SuboptCfg cfg);

}  // namespace mrna::brute

#endif  // COMPUTE_BRUTE_ALG_H_
