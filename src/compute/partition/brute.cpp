// Copyright 2021 Eliot Courtney.
#include "compute/partition/brute.h"

#include <utility>

#include "compute/brute.h"
#include "compute/partition/partition.h"
#include "model/primary.h"

namespace mrna::partition {

PartitionResult PartitionBruteForce(Primary r, const energy::EnergyModel& em) {
  // Allow lonely pairs for the partition function. TODO?
  auto res = BruteForce().Run(std::move(r), em, 1, true, false);
  return PartitionResult{.p = std::move(res.partition), .prob = std::move(res.probabilities)};
}

}  // namespace mrna::partition
