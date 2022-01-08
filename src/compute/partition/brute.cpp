// Copyright 2021 E.
#include "compute/partition/brute.h"

#include <set>
#include <stack>
#include <utility>
#include <vector>

#include "compute/brute.h"
#include "compute/energy/structure.h"
#include "compute/mfe/globals.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/partition.h"
#include "model/parsing.h"
#include "util/splaymap.h"

namespace mrna::partition {

std::pair<partition::Partition, partition::Probabilities> PartitionBruteForce(
    const Primary& r, const energy::EnergyModel& em) {
  // Allow lonely pairs for the partition function. TODO?
  auto res = BruteForce().Run(r, em, 1, true, false);
  return {std::move(res.partition), std::move(res.probabilities)};
}

}  // namespace mrna::partition
