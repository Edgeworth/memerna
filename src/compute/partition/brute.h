// Copyright 2021 Eliot Courtney.
#ifndef COMPUTE_PARTITION_BRUTE_H_
#define COMPUTE_PARTITION_BRUTE_H_

#include "compute/partition/partition.h"

namespace mrna::partition {

PartitionResult PartitionBruteForce(Primary r, const energy::EnergyModel& em);

}  // namespace mrna::partition

#endif  // COMPUTE_PARTITION_BRUTE_H_
