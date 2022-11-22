// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_BRUTE_ALG_H_
#define COMPUTE_BRUTE_ALG_H_

#include <vector>

#include "compute/energy/energy.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "compute/subopt/subopt_cfg.h"
#include "model/primary.h"

namespace mrna::brute {

subopt::SuboptResult MfeBruteForce(const Primary& r, energy::EnergyModelPtr em);

part::PartResult PartitionBruteForce(const Primary& r, energy::EnergyModelPtr em);

std::vector<subopt::SuboptResult> SuboptimalBruteForce(
    const Primary& r, energy::EnergyModelPtr em, subopt::SuboptCfg cfg);

}  // namespace mrna::brute

#endif  // COMPUTE_BRUTE_ALG_H_
