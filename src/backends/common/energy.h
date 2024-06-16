// Copyright 2024 Eliot Courtney.
#ifndef BACKENDS_COMMON_ENERGY_H_
#define BACKENDS_COMMON_ENERGY_H_

#include "model/energy.h"

namespace mrna::md {

inline Energy MinEnergy(const Energy* energy, std::size_t size) {
  Energy min = energy[0];
  for (int i = 0; i < static_cast<int>(size / sizeof(Energy)); ++i) min = std::min(min, energy[i]);
  return min;
}

}  // namespace mrna::md

#endif  // BACKENDS_COMMON_ENERGY_H_
