// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_MFE_BRUTE_H_
#define COMPUTE_MFE_BRUTE_H_

#include <utility>
#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"

namespace mrna::mfe {

Computed MfeBruteForce(const Primary& r, const energy::EnergyModel& em);

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_BRUTE_H_
