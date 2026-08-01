// Copyright 2024 Eliot Courtney.
#ifndef BACKENDS_COMMON_ENERGY_H_
#define BACKENDS_COMMON_ENERGY_H_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <type_traits>

#include "model/energy.h"
#include "util/util.h"

namespace mrna::md {

using ExpansionIndex = uint32_t;
constexpr ExpansionIndex INVALID_EXPANSION_INDEX = std::numeric_limits<ExpansionIndex>::max();

inline Energy MinEnergy(std::span<const Energy> energy) {
  return *std::min_element(energy.begin(), energy.end());
}

template <typename T>
  requires std::is_array_v<T>
inline Energy MinEnergy(const T& energy) {
  return MinEnergy(std::span<const Energy>(Decay(energy), sizeof(energy) / sizeof(Energy)));
}

}  // namespace mrna::md

#endif  // BACKENDS_COMMON_ENERGY_H_
