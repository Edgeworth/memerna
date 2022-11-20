// Copyright 2022 Eliot Courtney.
#include "model/energy.h"

#include <iomanip>

#include "model/constants.h"
#include "model/primary.h"
#include "util/error.h"

namespace mrna {

Energy Energy::FromDouble(double energy) {
  auto rounded = std::round(energy * Energy::FACTOR);
  auto res = Energy::FromRaw(static_cast<int32_t>(rounded));
  verify(res < CAP_E, "energy value out of range: %f", energy);
  return res;
}

Energy Energy::FromString(const std::string& s) {
  // Check for MAX and CAP.
  if (s == "MAX") return MAX_E;
  if (s == "CAP") return CAP_E;
  // Check for a number.
  std::stringstream ss(s);
  int sgn = ss.peek() == '-' ? -1 : 1;
  int hi = 0;
  int lo = 0;
  verify(ss >> hi, "invalid energy: %s", s.c_str());
  hi = std::abs(hi);
  if (!ss.eof()) {
    verify(ss.get() == '.', "invalid energy: %s", s.c_str());
    std::string decimal;
    verify(ss >> decimal, "invalid energy: %s", s.c_str());
    verify(decimal.size() <= Energy::EXPONENT, "invalid energy: %s", s.c_str());
    decimal.resize(Energy::EXPONENT, '0');
    lo = std::stoi(decimal);
  }
  verify(ss.eof(), "invalid energy: %s", s.c_str());
  verify(lo >= 0 && lo < Energy::FACTOR, "invalid energy: %s", s.c_str());
  Energy energy = Energy::FromRaw(sgn * (hi * Energy::FACTOR + lo));
  verify(energy < CAP_E, "energy out of range: %s", s.c_str());
  return energy;
}

std::string Energy::ToString() const {
  if (*this == MAX_E) return "MAX";
  if (*this == CAP_E) return "CAP";
  std::stringstream ss;
  const auto* sgn = (v < 0 ? "-" : "");
  const auto abs = std::abs(v);
  ss << sgn << abs / FACTOR << "." << std::setw(EXPONENT) << std::setfill('0') << abs % FACTOR;
  return ss.str();
}

BoltzEnergy Energy::Boltz() const {
  if (*this >= CAP_E) return 0;
  return exp(BoltzEnergy(ToDouble()) * (BoltzEnergy(-1) / BoltzEnergy(R * T)));
}

std::istream& operator>>(std::istream& str, Energy& o) {
  std::string s;
  str >> s;
  o = Energy::FromString(s);
  return str;
}

std::ostream& operator<<(std::ostream& out, const Energy& o) {
  out << o.ToString();
  return out;
}

}  // namespace mrna
