// Copyright 2022 Eliot Courtney.
#include "model/model.h"

#include <iomanip>

#include "model/primary.h"
#include "util/error.h"

namespace mrna {

Energy EnergyFromString(const std::string& s) {
  // Check for MAX and CAP.
  if (s == "MAX") return MAX_E;
  if (s == "CAP") return CAP_E;
  // Check for a number.
  std::stringstream ss(s);
  int hi = 0;
  int lo = 0;
  verify(ss >> hi, "invalid energy: %s", s.c_str());
  if (!ss.eof()) {
    verify(ss.get() == '.', "invalid energy: %s", s.c_str());
    verify(ss >> lo, "invalid energy: %s", s.c_str());
  }
  verify(ss.eof(), "invalid energy: %s", s.c_str());
  verify(lo >= 0 && lo < ENERGY_FACTOR, "invalid energy: %s", s.c_str());
  Energy energy = hi * ENERGY_FACTOR + lo;
  verify(energy < CAP_E && energy > -CAP_E, "energy out of range: %s", s.c_str());
  return energy;
}

std::string EnergyToString(Energy energy) {
  if (energy == MAX_E) return "MAX";
  if (energy == CAP_E) return "CAP";
  std::stringstream ss;
  ss << energy / ENERGY_FACTOR << "." << std::setw(ENERGY_EXPONENT) << std::setfill('0')
     << abs(energy % ENERGY_FACTOR);
  return ss.str();
}

}  // namespace mrna
