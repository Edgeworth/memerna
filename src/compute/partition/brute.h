// Copyright 2021 Eliot Courtney.
#ifndef COMPUTE_PARTITION_BRUTE_H_
#define COMPUTE_PARTITION_BRUTE_H_

#include <utility>
#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"

namespace mrna::partition {

std::pair<partition::partition_t, partition::probabilities_t> PartitionBruteForce(
    const primary_t& r, const energy::EnergyModel& em);

}  // namespace mrna::partition

#endif  // COMPUTE_PARTITION_BRUTE_H_
