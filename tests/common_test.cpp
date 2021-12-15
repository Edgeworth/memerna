// Copyright 2016 E.
#include "common_test.h"

#include "energy/structure.h"
#include "model/parsing.h"

namespace mrna {

energy::EnergyModelPtr g_em;

std::ostream& operator<<(std::ostream& os, const secondary_t& s) {
  return os << "(" << PrimaryToString(s.r) << ", " << PairsToDotBracket(s.p) << ")";
}

std::ostream& operator<<(std::ostream& os, const computed_t& computed) {
  os << computed.s;
  os << ", (";
  for (auto ctd : computed.base_ctds) os << energy::CtdToName(ctd) << ", ";
  os << "), " << computed.energy;
  return os;
}

}  // namespace mrna
