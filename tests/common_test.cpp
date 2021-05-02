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
#include "common_test.h"

#include "energy/structure.h"
#include "parsing.h"

namespace memerna {

energy::EnergyModelPtr g_em;
std::vector<energy::EnergyModelPtr> g_ems;

std::ostream& operator<<(std::ostream& os, const secondary_t& s) {
  return os << "(" << parsing::PrimaryToString(s.r) << ", " << parsing::PairsToDotBracket(s.p)
            << ")";
}

std::ostream& operator<<(std::ostream& os, const computed_t& computed) {
  os << computed.s;
  os << ", (";
  for (auto ctd : computed.base_ctds) os << energy::CtdToName(ctd) << ", ";
  os << "), " << computed.energy;
  return os;
}
}  // namespace memerna
