// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_INTERNAL_H_
#define COMPUTE_ENERGY_INTERNAL_H_

#include <deque>
#include <utility>

#include "common.h"
#include "compute/energy/model.h"

namespace mrna::energy::internal {

typedef std::deque<std::pair<Ctd, energy_t>> branch_ctd_t;

energy_t ComputeOptimalCtd(const secondary_t& secondary, const EnergyModel& em,
    const std::deque<int>& branches, bool use_first_lu, branch_ctd_t& branch_ctds);

// Per-base representation:
// At the index of the right side of a branch, we store the CTD identifier
// for this branch for when it is an outer loop (i.e. closing a multiloop).
// At the index of the left side, we store the CTD identifier for this branch for when it is inside
// a multiloop. This is because in this situation ((...).(...))(...), the first loop can coaxially
// stack twice, one when it is an outer loop, and the other with the loop on the far right.

// List representation:
// By convention, if there is an outer loop, it comes first (not last), but it does not matter to
// these functions.
// Outer loops are represented by branch indices referring to the right side of the loop (i.e.
// p[branch] < branch).

// Takes the list representation of ctds in |branch_ctds| for |branches| branches and
// writes it in per-base representation to |computed|.
void AddBranchCtdsToComputed(
    computed_t& computed, const std::deque<int>& branches, const branch_ctd_t& branch_ctds);
// Reads the per-base ctd representation from |computed| for |branches| branches and
// writes it in list representation to |branch_ctds|.
energy_t GetBranchCtdsFromComputed(const computed_t& computed, const EnergyModel& em,
    const std::deque<int>& branches, branch_ctd_t& branch_ctds);

}  // namespace mrna::energy::internal

#endif  // COMPUTE_ENERGY_INTERNAL_H_
