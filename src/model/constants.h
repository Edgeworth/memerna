// Copyright 2022 Eliot Courtney.
#ifndef MODEL_CONSTANTS_H_
#define MODEL_CONSTANTS_H_

#include "model/energy.h"

namespace mrna {

// -----------------------------------------------
// Values affecting the energy model:
inline constexpr int HAIRPIN_MIN_SZ = 3;
// N.B. This is for kcal/mol so it's not 8.315.
inline constexpr double R = 1.98720425864083e-3;
// This is 37 degrees Celsius. Changing this is not a good idea.
inline constexpr double T = 310.15;
// Ninio maximum asymmetry.
inline const Energy NINIO_MAX_ASYM = E(3.0);
// Maximum size of a twoloop.
inline constexpr int TWOLOOP_MAX_SZ = 30;

}  // namespace mrna

#endif  // MODEL_CONSTANTS_H_
