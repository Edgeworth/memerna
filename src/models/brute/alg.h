// Copyright 2022 Eliot Courtney.
#ifndef MODELS_BRUTE_ALG_H_
#define MODELS_BRUTE_ALG_H_

#include <vector>

#include "api/energy/model.h"
#include "api/part.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "model/primary.h"

namespace mrna::md::brute {

subopt::SuboptResult MfeBrute(const Primary& r, erg::EnergyModelPtr em);

part::PartResult PartitionBrute(const Primary& r, erg::EnergyModelPtr em);

std::vector<subopt::SuboptResult> SuboptBrute(
    const Primary& r, erg::EnergyModelPtr em, subopt::SuboptCfg cfg);

}  // namespace mrna::md::brute

#endif  // MODELS_BRUTE_ALG_H_
