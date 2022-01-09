// Copyright 2016 Eliot Courtney.
#include "common_test.h"

#include "compute/energy/structure.h"

namespace mrna {

energy::EnergyModel g_em;

std::ostream& operator<<(std::ostream& os, const Computed& c) {
  os << "(" << PrimaryToString(c.r) << ", " << SecondaryToDotBracket(c.s) << ")";
  os << ", (";
  for (auto ctd : c.base_ctds) os << energy::CtdToName(ctd) << ", ";
  os << "), " << c.energy;
  return os;
}

}  // namespace mrna
