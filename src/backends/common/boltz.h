// Copyright 2023 Eliot Courtney.
#ifndef BACKENDS_COMMON_BOLTZ_H_
#define BACKENDS_COMMON_BOLTZ_H_

#include "model/energy.h"
#include "util/util.h"  // IWYU pragma: keep

namespace mrna::md {

inline void FillBoltzArray(BoltzEnergy* output, const Energy* input, int elements) {
  for (int i = 0; i < elements; ++i) output[i] = input[i].Boltz();
}

#define FILL_BOLTZ(bm, name)                                                                  \
  /* NOLINTBEGIN */                                                                           \
  static_assert(sizeof((bm).name) / sizeof(*Decay((bm).name)) ==                              \
          sizeof((bm).m().name) / sizeof(*Decay((bm).m().name)),                              \
      "BoltzModel does not match Model");                                                     \
  FillBoltzArray(                                                                             \
      Decay((bm).name), Decay((bm).m().name), sizeof((bm).name) / sizeof(*Decay((bm).name))); \
  /* NOLINTEND */

}  // namespace mrna::md

#endif  // BACKENDS_COMMON_BOLTZ_H_
