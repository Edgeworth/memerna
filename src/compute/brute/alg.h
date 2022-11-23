// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_BRUTE_ALG_H_
#define COMPUTE_BRUTE_ALG_H_

#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "compute/subopt/subopt_cfg.h"
#include "model/primary.h"

namespace mrna::brute {

subopt::SuboptResult MfeBrute(const Primary& r, erg::EnergyModelPtr em);

part::PartResult PartitionBrute(const Primary& r, erg::EnergyModelPtr em);

std::vector<subopt::SuboptResult> SuboptBrute(
    const Primary& r, erg::EnergyModelPtr em, subopt::SuboptCfg cfg);

}  // namespace mrna::brute

#endif  // COMPUTE_BRUTE_ALG_H_
