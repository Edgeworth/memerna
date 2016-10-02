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
