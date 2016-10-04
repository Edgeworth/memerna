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

namespace memerna {
namespace fold {
namespace internal {

primary_t gr;
std::vector<int> gp;
std::vector<Ctd> gctd;
std::string grep;
precomp_t gpc;
energy_t genergy;
energy::EnergyModel gem;
array3d_t<energy_t, DP_SIZE> gdp;
array2d_t<energy_t, EXT_SIZE> gext;

void SetGlobalState(const primary_t& r, const energy::EnergyModel& em) {
  gr = r;
  gp.resize(gr.size());
  gctd.resize(gr.size());
  genergy = MAX_E;
  gem = em;
  gpc = PrecomputeData(gr, gem);
  std::fill(gp.begin(), gp.end(), -1);
  std::fill(gctd.begin(), gctd.end(), CTD_NA);
  gdp = array3d_t<energy_t, DP_SIZE>(gr.size() + 1);
  gext = array2d_t<energy_t, EXT_SIZE>(gr.size() + 1);
}
}
}
}
