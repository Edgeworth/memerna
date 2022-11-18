// Copyright 2022 Eliot Courtney.
#ifndef MODEL_ENERGY_H_
#define MODEL_ENERGY_H_

#include <algorithm>
#include <string>

#include "model/ctd.h"
#include "model/primary.h"
#include "util/error.h"
#include "util/float.h"

namespace mrna {

using BoltzEnergy = flt;

// Make sure it is packed and 4 byte aligned so memset works.
struct __attribute__((packed, aligned(4))) Energy {
 public:
  // Precision of energy values.
  static constexpr int FACTOR = 100;
  static constexpr int EXPONENT = 2;

  [[nodiscard]] static constexpr Energy FromRaw(int32_t v) { return Energy{.v = v}; }

  // Converts a floating point energy value in kcal/mol to an integer energy value.
  [[nodiscard]] static Energy FromDouble(double energy);

  [[nodiscard]] static Energy FromString(const std::string& s);

  [[nodiscard]] std::string ToString() const;

  [[nodiscard]] double ToDouble() const { return v / static_cast<double>(FACTOR); }

  [[nodiscard]] BoltzEnergy Boltz() const;

  constexpr auto operator<=>(const Energy&) const = default;

  constexpr Energy operator-() const { return FromRaw(-v); }

  constexpr Energy operator+(const Energy& o) const { return FromRaw(v + o.v); }
  constexpr Energy operator+=(const Energy& o) {
    v += o.v;
    return *this;
  }

  constexpr Energy operator-(const Energy& o) const { return FromRaw(v - o.v); }
  constexpr Energy operator-=(const Energy& o) {
    v -= o.v;
    return *this;
  }

  int32_t v;
};

constexpr Energy operator*(const Energy& e, int o) { return Energy::FromRaw(e.v * o); }
constexpr Energy operator*(int o, const Energy& e) { return Energy::FromRaw(e.v * o); }

std::istream& operator>>(std::istream& str, Energy& o);
std::ostream& operator<<(std::ostream& out, const Energy& o);

// Don't change these values. Plays nice with memset.
// Used for infinite/sentinel energy values, e.g. in DP tables.
inline constexpr Energy MAX_E = Energy::FromRaw(0x0F0F0F0F);

// Used for finite but larger than any possible energy values. e.g. for subopt-delta
inline constexpr Energy CAP_E = Energy::FromRaw(0x07070707);

inline constexpr Energy ZERO_E = Energy::FromRaw(0);

[[nodiscard]] inline Energy E(double energy) { return Energy::FromDouble(energy); }

}  // namespace mrna

#endif  // MODEL_ENERGY_H_
