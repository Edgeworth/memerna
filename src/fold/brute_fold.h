// Copyright 2016 E.
#ifndef FOLD_BRUTE_FOLD_H_
#define FOLD_BRUTE_FOLD_H_

#include <utility>
#include <vector>

#include "common.h"
#include "energy/energy_model.h"
#include "partition/partition.h"

namespace mrna::fold {

namespace internal {
std::vector<int> GetBranchCounts(const std::vector<int>& p);
}

computed_t FoldBruteForce(const primary_t& r, const energy::EnergyModel& em);

std::vector<computed_t> SuboptimalBruteForce(
    const primary_t& r, const energy::EnergyModel& em, int max_structures_);

std::pair<partition::partition_t, partition::probabilities_t> PartitionBruteForce(
    const primary_t& r, const energy::EnergyModel& em);

}  // namespace mrna::fold
#endif  // FOLD_BRUTE_FOLD_H_
