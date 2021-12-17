// Copyright 2016 E.
#ifndef COMPUTE_MFE_BRUTE_H_
#define COMPUTE_MFE_BRUTE_H_

#include <utility>
#include <vector>

#include "common.h"
#include "compute/energy/model.h"
#include "compute/partition/partition.h"

namespace mrna::mfe {

computed_t MfeBruteForce(const primary_t& r, const energy::EnergyModel& em);

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_BRUTE_H_
