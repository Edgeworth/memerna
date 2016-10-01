#ifndef MEMERNA_BRUTE_FOLD_H
#define MEMERNA_BRUTE_FOLD_H

#include "common.h"
#include "energy/energy_model.h"

namespace memerna {
namespace fold {

namespace internal {
std::vector<int> GetBranchCounts(const std::vector<int>& p);
}

std::vector<computed_t> FoldBruteForce(
    const primary_t& r, const energy::EnergyModel& em, int max_structures_);
}
}
#endif  // MEMERNA_BRUTE_FOLD_H
