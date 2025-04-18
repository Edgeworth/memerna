// Copyright 2022 E.
#ifndef MODEL_CONSTANTS_H_
#define MODEL_CONSTANTS_H_

#include "util/float.h"

namespace mrna {

// -----------------------------------------------
// Values affecting the energy model:
inline constexpr int HAIRPIN_MIN_SZ = 3;
// N.B. This is for kcal/mol so it's not 8.315.
inline const flt R = FLT(1.98720425864083e-3);
// This is 37 degrees Celsius. Changing this is not a good idea.
inline const flt T = FLT(310.15);
// Maximum size of a twoloop.
inline constexpr int TWOLOOP_MAX_SZ = 30;

}  // namespace mrna

#endif  // MODEL_CONSTANTS_H_
