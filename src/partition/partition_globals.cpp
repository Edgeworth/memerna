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
#include "partition/partition_globals.h"
#include "energy/energy_globals.h"

namespace memerna {
namespace partition {
namespace internal {

array3d_t<penergy_t, PT_SIZE> gpt;
array2d_t<penergy_t, PTEXT_SIZE> gptext;
precomp_t gppc;

}

void SetPartitionGlobalState(const primary_t& r, const energy::EnergyModel& em) {
  energy::SetEnergyGlobalState(r, em);
  // 0.0 is zero'd memory.
  internal::gpt = array3d_t<penergy_t, PT_SIZE>(gr.size() + 1, 0);
  internal::gptext = array2d_t<penergy_t, PTEXT_SIZE>(gr.size() + 1, 0);
  internal::gppc = internal::PrecomputeData(r, em);
}

}
}
