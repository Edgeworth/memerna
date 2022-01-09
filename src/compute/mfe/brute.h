// Copyright 2016 E.
#ifndef COMPUTE_MFE_BRUTE_H_
#define COMPUTE_MFE_BRUTE_H_

#include <utility>
#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"

namespace mrna::mfe {

subopt::SuboptResult MfeBruteForce(const Primary& r, const energy::EnergyModel& em);

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_BRUTE_H_
