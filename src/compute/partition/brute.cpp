// Copyright 2021 E.
#include "compute/partition/brute.h"

#include <set>
#include <stack>
#include <utility>
#include <vector>

#include "compute/brute.h"
#include "compute/energy/structure.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/partition.h"
#include "util/splaymap.h"

namespace mrna::partition {

PartitionResult PartitionBruteForce(const Primary& r, const energy::EnergyModel& em) {
  // Allow lonely pairs for the partition function. TODO?
  auto res = BruteForce().Run(r, em, 1, true, false);
  return PartitionResult{.p = std::move(res.partition), .prob = std::move(res.probabilities)};
}

}  // namespace mrna::partition
