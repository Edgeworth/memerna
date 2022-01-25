// Copyright 2016 E.
#ifndef COMPUTE_MFE_BRUTE_H_
#define COMPUTE_MFE_BRUTE_H_

#include "compute/energy/model.h"
#include "compute/subopt/subopt.h"

namespace mrna::mfe {

subopt::SuboptResult MfeBruteForce(Primary r, const energy::EnergyModel& em);

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_BRUTE_H_
