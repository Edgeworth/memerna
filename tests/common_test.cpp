// Copyright 2016 Eliot Courtney.
#include "common_test.h"

#include "compute/energy/structure.h"
#include "model/parsing.h"

namespace mrna {

energy::EnergyModel g_em;

std::ostream& operator<<(std::ostream& os, const Secondary& s) {
  return os << "(" << PrimaryToString(s.r) << ", " << PairsToDotBracket(s.p) << ")";
}

std::ostream& operator<<(std::ostream& os, const Computed& computed) {
  os << computed.s;
  os << ", (";
  for (auto ctd : computed.base_ctds) os << energy::CtdToName(ctd) << ", ";
  os << "), " << computed.energy;
  return os;
}

}  // namespace mrna
