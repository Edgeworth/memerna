// Copyright 2022 Eliot Courtney.
#include "model/energy.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <istream>
#include <random>

#include "util/error.h"

namespace mrna {

Energy Energy::FromFlt(flt energy) {
  auto rounded = round(energy * Energy::FACTOR);
  auto res = Energy::FromRaw(static_cast<int32_t>(rounded));
  verify(res < CAP_E, "energy value out of range: {}", energy);
  return res;
}

Energy Energy::FromString(const std::string& s) {
  // Check for MAX and CAP.
  if (s == "MAX") return MAX_E;
  if (s == "CAP") return CAP_E;
  // Check for a number.
  std::stringstream ss(s);
  const int sgn = ss.peek() == '-' ? -1 : 1;
  int hi = 0;
  int lo = 0;
  verify(ss >> hi, "invalid energy: {}", s);
  hi = std::abs(hi);
  if (!ss.eof()) {
    verify(ss.get() == '.', "invalid energy: {}", s);
    std::string decimal;
    verify(ss >> decimal, "invalid energy: {}", s);
    verify(decimal.size() <= Energy::EXPONENT, "invalid energy: {}", s);
    decimal.resize(Energy::EXPONENT, '0');
    lo = std::stoi(decimal);
  }
  verify(ss.eof(), "invalid energy: {}", s);
  verify(lo >= 0 && lo < Energy::FACTOR, "invalid energy: {}", s);
  Energy energy = Energy::FromRaw(sgn * (hi * Energy::FACTOR + lo));
  verify(energy < CAP_E, "energy out of range: {}", s);
  return energy;
}

std::string Energy::ToString() const noexcept {
  if (*this == MAX_E) return "MAX";
  if (*this == CAP_E) return "CAP";
  std::stringstream ss;
  const auto* sgn = (v < 0 ? "-" : "");
  const auto abs = std::abs(v);
  ss << sgn << abs / FACTOR << "." << std::setw(EXPONENT) << std::setfill('0') << abs % FACTOR;
  return ss.str();
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

std::vector<Energy> RandomEnergies(
    std::size_t length, Energy min_energy, Energy max_energy, std::mt19937& eng) {
  std::uniform_int_distribution<decltype(min_energy.v)> energy_dist(min_energy.v, max_energy.v);
  std::vector<Energy> energies(length);
  for (int i = 0; i < int(length); ++i) energies[i] = Energy::FromRaw(energy_dist(eng));
  return energies;
}

}  // namespace mrna
