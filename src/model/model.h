// Copyright 2022 Eliot Courtney.
#ifndef MODEL_MODEL_H_
#define MODEL_MODEL_H_

#include <algorithm>

#include "model/primary.h"
#include "util/error.h"
#include "util/float.h"

namespace mrna {

using Energy = int32_t;
using BoltzEnergy = flt;

// Don't change these values. Plays nice with memset.
// Used for infinite/sentinel energy values, e.g. in DP tables.
inline constexpr Energy MAX_E = 0x0F0F0F0F;
// Used for finite but larger than any possible energy values. e.g. for subopt-delta
inline constexpr Energy CAP_E = 0x07070707;

// Precision of energy values.
inline constexpr int ENERGY_FACTOR = 100;
inline constexpr int ENERGY_EXPONENT = 2;

// Converts a floating point energy value in kcal/mol to an integer energy value.
inline Energy E(double energy) {
  auto rounded = static_cast<Energy>(std::round(energy * ENERGY_FACTOR));
  verify(abs_eq(energy * ENERGY_FACTOR, rounded), "energy value not convertible from double: %f",
      energy);
  verify(rounded < CAP_E && rounded > -CAP_E, "energy value out of range: %f", energy);
  return rounded;
}

Energy EnergyFromString(const std::string& s);

std::string EnergyToString(Energy energy);

// -----------------------------------------------
// Values affecting the energy model:
inline constexpr int HAIRPIN_MIN_SZ = 3;
// N.B. This is for kcal/mol so it's not 8.315.
inline constexpr double R = 1.98720425864083e-3;
// This is 37 degrees Celsius. Changing this is not a good idea.
inline constexpr double T = 310.15;
// Ninio maximum asymmetry.
inline constexpr Energy NINIO_MAX_ASYM = 30;
// Maximum size of a twoloop.
inline constexpr int TWOLOOP_MAX_SZ = 30;

inline BoltzEnergy Boltz(Energy energy) {
  if (energy >= CAP_E) return 0;
  return exp(BoltzEnergy(energy) * (BoltzEnergy(-1) / BoltzEnergy(ENERGY_FACTOR * R * T)));
}

}  // namespace mrna

#endif  // MODEL_MODEL_H_
