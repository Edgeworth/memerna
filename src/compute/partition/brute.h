// Copyright 2021 E.
#ifndef COMPUTE_PARTITION_BRUTE_H_
#define COMPUTE_PARTITION_BRUTE_H_

#include <utility>
#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"

namespace mrna::partition {

std::pair<partition::Partition, partition::Probabilities> PartitionBruteForce(
    const Primary& r, const energy::EnergyModel& em);

}  // namespace mrna::partition

#endif  // COMPUTE_PARTITION_BRUTE_H_
