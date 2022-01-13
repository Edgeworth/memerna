// Copyright 2016 E.
#ifndef COMPUTE_ENERGY_INTERNAL_H_
#define COMPUTE_ENERGY_INTERNAL_H_

#include <deque>
#include <utility>

#include "compute/energy/model.h"
#include "model/ctd.h"

namespace mrna::energy::internal {

using BranchCtd = std::deque<std::pair<Ctd, Energy>>;

Energy ComputeOptimalCtds(const Primary& r, const Secondary& s, const EnergyModel& em,
    const std::deque<int>& branches, bool use_first_lu, BranchCtd* branch_ctd);

// Per-base representation of CTDs:
// An array the same length as the sequence, which contains CTD identifiers.
// For each base-pair at a branch, we store a CTD identifier or nothing.
// Since there are two locations per base-pair, we can store up to two things.
// At the index of the right side of a branch, we store the CTD identifier
// for this branch for when it is an outer loop (i.e. closing a multiloop).
// At the index of the left side, we store the CTD identifier for this branch
// for when it is inside a multiloop. This is because in this situation
// ((...).(...))(...), the first loop can coaxially stack twice, one when it is
// an outer loop, and the other with the loop on the far right.

// Branch representation:
// A list of indices of branches. Depending on which side of the branch's index
// is used, it refers to the branch as an outer loop or an inner loop.
// By convention, if there is an outer loop, it comes first (not last), but it
// does not matter to these functions. Outer loops are represented by branch
// indices referring to the right side of the loop (i.e. s[branch] < branch).

// Takes the branch representation of ctds in |branch_ctd| for |branches| branches and
// writes it in per-base representation to |ctds|.
void AddBranchCtdsToBaseCtds(
    const std::deque<int>& branches, const BranchCtd& branch_ctd, Ctds* ctd);

// Reads the per-base ctd representation from |ctd| for |branches| branches and
// writes it in branch representation to |branch_ctd|.
Energy AddBaseCtdsToBranchCtds(const Primary& r, const Secondary& s, const Ctds& ctd,
    const EnergyModel& em, const std::deque<int>& branches, BranchCtd* branch_ctd);

}  // namespace mrna::energy::internal

#endif  // COMPUTE_ENERGY_INTERNAL_H_
