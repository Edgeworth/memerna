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
#include "globals.h"
#include "fold/fold_globals.h"
#include "energy/energy_globals.h"

namespace memerna {
namespace fold {
namespace internal {

std::vector<int> gp;
std::vector<Ctd> gctd;
std::string grep;
energy_t genergy;
array3d_t<energy_t, DP_SIZE> gdp;
array2d_t<energy_t, EXT_SIZE> gext;

}

void SetFoldGlobalState(const primary_t& r, const energy::EnergyModel& em) {
  energy::SetEnergyGlobalState(r, em);
  internal::gp.resize(gr.size());
  internal::gctd.resize(gr.size());
  internal::genergy = MAX_E;
  std::fill(internal::gp.begin(), internal::gp.end(), -1);
  std::fill(internal::gctd.begin(), internal::gctd.end(), CTD_NA);
  internal::gdp = array3d_t<energy_t, internal::DP_SIZE>(gr.size() + 1);
  internal::gext = array2d_t<energy_t, internal::EXT_SIZE>(gr.size() + 1);
}

}
}
