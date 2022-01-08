#ifndef MODEL_MODEL_H_
#define MODEL_MODEL_H_

#include <cstdint>

#include "model/primary.h"
#include "util/float.h"

namespace mrna {

typedef int32_t Energy;
typedef flt BoltzEnergy;

// Don't change this value. Plays nice with memset.
inline constexpr Energy MAX_E = 0x0F0F0F0F;
inline constexpr Energy CAP_E = 0x07070707;

// -----------------------------------------------
// Values affecting the energy model:
inline constexpr int HAIRPIN_MIN_SZ = 3;
// N.B. This is for kcal/mol so it's not 8.315.
inline constexpr double R = 1.9872036e-3;
// This is 37 degrees Celsius. Changing this is not a good idea.
inline constexpr double T = 310.15;
// Ninio maximum asymmetry.
inline constexpr Energy NINIO_MAX_ASYM = 30;
// Maximum size of a twoloop.
inline constexpr int TWOLOOP_MAX_SZ = 30;

inline bool ViableFoldingPair(const Primary& r, int st, int en) {
  return CanPair(r[st], r[en]) && (en - st - 1 >= HAIRPIN_MIN_SZ) &&
      ((en - st - 3 >= HAIRPIN_MIN_SZ && CanPair(r[st + 1], r[en - 1])) ||
          (st > 0 && en < static_cast<int>(r.size() - 1) && CanPair(r[st - 1], r[en + 1])));
}

inline BoltzEnergy Boltz(Energy energy) {
  if (energy >= CAP_E) return 0;
  return exp(BoltzEnergy(energy) * (BoltzEnergy(-1) / BoltzEnergy(10.0 * R * T)));
}

}  // namespace mrna

#endif  // MODEL_MODEL_H_
