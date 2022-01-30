// Copyright 2021 E.
#include "compute/partition/brute.h"

#include <utility>

#include "compute/brute/brute.h"
#include "compute/partition/partition.h"
#include "model/primary.h"

namespace mrna::partition {

PartitionResult PartitionBruteForce(Primary r, const energy::EnergyModel& em) {
  // TODO: Allow lonely pairs for the partition function
  auto res = brute::BruteForce().Run(std::move(r), em, 1, true, false);
  return PartitionResult{
      .dp{}, .ext{}, .p = std::move(res.partition), .prob = std::move(res.probabilities)};
}

}  // namespace mrna::partition
