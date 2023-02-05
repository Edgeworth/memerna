// Copyright 2023 Eliot Courtney.
#ifndef COMPUTE_ENERGY_COMMON_BOLTZ_H_
#define COMPUTE_ENERGY_COMMON_BOLTZ_H_

#include "model/energy.h"
#include "util/util.h"

namespace mrna::md {

inline void FillBoltzArray(BoltzEnergy* output, const Energy* input, int elements) {
  for (int i = 0; i < elements; ++i) output[i] = input[i].Boltz();
}

#define FILL_BOLTZ(name)                                                                    \
  static_assert(/* NOLINTNEXTLINE */                                                        \
      sizeof(name) / sizeof(*Decay(name)) == sizeof(em().name) / sizeof(*Decay(em().name)), \
      "BoltzModel does not match Model");                                                   \
  /* NOLINTNEXTLINE */                                                                      \
  FillBoltzArray(Decay(name), Decay(em().name), sizeof(name) / sizeof(*Decay(name)));

}  // namespace mrna::md

#endif  // COMPUTE_ENERGY_COMMON_BOLTZ_H_
