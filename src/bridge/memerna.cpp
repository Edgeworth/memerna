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
#include "bridge/memerna.h"
#include "energy/structure.h"

namespace memerna {
namespace bridge {

energy_t Memerna::Efn(const secondary_t& secondary, std::string* desc) const {
  computed_t computed;
  if (desc) {
    std::unique_ptr<energy::Structure> structure;
    computed = energy::ComputeEnergy(secondary, *em, &structure);
    for (const auto& s : structure->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    computed = energy::ComputeEnergy(secondary, *em);
  }

  return computed.energy;
}

computed_t Memerna::Fold(const primary_t& r) const { return fold::Context(r, em, options).Fold(); }

int Memerna::Suboptimal(fold::SuboptimalCallback fn,
    const primary_t& r, energy_t energy_delta) const {
  return fold::Context(r, em, options).Suboptimal(fn, true, energy_delta, -1);
}

std::vector<computed_t>
Memerna::SuboptimalIntoVector(const primary_t& r, energy_t energy_delta) const {
  return fold::Context(r, em, options).SuboptimalIntoVector(true, energy_delta, -1);
}
}
}
