// Copyright 2023 Eliot Courtney.
#ifndef MODELS_COMMON_BOLTZ_H_
#define MODELS_COMMON_BOLTZ_H_

#include "model/energy.h"
#include "util/util.h"

namespace mrna::md {

inline void FillBoltzArray(BoltzEnergy* output, const Energy* input, int elements) {
  for (int i = 0; i < elements; ++i) output[i] = input[i].Boltz();
}

#define FILL_BOLTZ(bem, name)                                        \
  static_assert(/* NOLINTNEXTLINE */                                 \
      sizeof(bem.name) / sizeof(*Decay(bem.name)) ==                 \
          sizeof((bem).em().name) / sizeof(*Decay((bem).em().name)), \
      "BoltzModel does not match Model");                            \
  /* NOLINTNEXTLINE */                                               \
  FillBoltzArray(                                                    \
      Decay((bem).name), Decay((bem).em().name), sizeof((bem).name) / sizeof(*Decay((bem).name)));

}  // namespace mrna::md

#endif  // MODELS_COMMON_BOLTZ_H_
