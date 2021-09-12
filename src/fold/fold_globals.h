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
#ifndef MEMERNA_FOLD_GLOBALS_H_
#define MEMERNA_FOLD_GLOBALS_H_

#include "array.h"
#include "common.h"
#include "energy/energy_model.h"
#include "fold/fold_constants.h"

namespace memerna {
namespace fold {
namespace internal {

extern std::vector<int> gp;
extern std::vector<Ctd> gctd;
extern std::string grep;
extern energy_t genergy;
extern array3d_t<energy_t, DP_SIZE> gdp;
extern array2d_t<energy_t, EXT_SIZE> gext;

}  // namespace internal

void SetFoldGlobalState(const primary_t& r, const energy::EnergyModel& em);

}  // namespace fold
}  // namespace memerna

#endif  // MEMERNA_FOLD_GLOBALS_H_
