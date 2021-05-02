// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#ifndef MEMERNA_ENERGY_INTERNAL_H
#define MEMERNA_ENERGY_INTERNAL_H

#include <deque>

#include "common.h"
#include "energy/energy_model.h"

namespace memerna {
namespace energy {
namespace internal {

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
}  // namespace internal
}  // namespace energy
}  // namespace memerna

#endif  // MEMERNA_ENERGY_INTERNAL_H
