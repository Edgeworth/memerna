// Copyright 2022 Eliot Courtney.
#ifndef MODEL_ENERGY_H_
#define MODEL_ENERGY_H_

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <compare>
#include <cstdint>
#include <iosfwd>
#include <string>

#include "model/constants.h"
#include "util/float.h"

namespace mrna {

using BoltzEnergy = flt;

// Make sure it is packed and 4 byte aligned so memset works.
struct __attribute__((packed, aligned(4))) Energy {
 public:
  // Precision of energy values.
  static constexpr int FACTOR = powi(10, ENERGY_PRECISION);
  static constexpr int EXPONENT = ENERGY_PRECISION;

  [[nodiscard]] inline static constexpr Energy FromRaw(int32_t v) noexcept {
    return Energy{.v = v};
  }

  // Converts a floating point energy value in kcal/mol to an integer energy value.
  [[nodiscard]] static Energy FromDouble(double energy);
  [[nodiscard]] static Energy FromString(const std::string& s);

  [[nodiscard]] std::string ToString() const noexcept;
  [[nodiscard]] inline double ToDouble() const noexcept { return v / static_cast<double>(FACTOR); }
  [[nodiscard]] inline BoltzEnergy Boltz() const noexcept;

  inline constexpr auto operator<=>(const Energy&) const noexcept = default;

  inline constexpr Energy operator-() const noexcept { return FromRaw(-v); }

  inline constexpr Energy operator+(const Energy& o) const noexcept { return FromRaw(v + o.v); }
  inline constexpr Energy operator+=(const Energy& o) noexcept {
    v += o.v;
    return *this;
  }

  inline constexpr Energy operator-(const Energy& o) const noexcept { return FromRaw(v - o.v); }
  inline constexpr Energy operator-=(const Energy& o) {
    v -= o.v;
    return *this;
  }

  int32_t v;
};

constexpr Energy operator*(const Energy& e, int o) noexcept { return Energy::FromRaw(e.v * o); }
constexpr Energy operator*(int o, const Energy& e) noexcept { return Energy::FromRaw(e.v * o); }

std::istream& operator>>(std::istream& str, Energy& o);
std::ostream& operator<<(std::ostream& out, const Energy& o);

[[nodiscard]] inline Energy E(double energy) { return Energy::FromDouble(energy); }

// Don't change these values. Plays nice with memset.
// Used for infinite/sentinel energy values, e.g. in DP tables.
inline constexpr Energy MAX_E = Energy::FromRaw(0x0F0F0F0F);

// Used for finite but larger than any possible energy values. e.g. for subopt-delta
inline constexpr Energy CAP_E = Energy::FromRaw(0x07070707);

inline constexpr Energy ZERO_E = Energy::FromRaw(0);

// Ninio maximum asymmetry.
inline const Energy NINIO_MAX_ASYM = E(3.0);

[[nodiscard]] inline BoltzEnergy Energy::Boltz() const noexcept {
  if (*this >= CAP_E) return 0;
  return exp(BoltzEnergy(-ToDouble()) / (BoltzEnergy(R) * BoltzEnergy(T)));
}

}  // namespace mrna

template <>
struct fmt::formatter<mrna::Energy> : ostream_formatter {};

#endif  // MODEL_ENERGY_H_
