// Copyright 2016 Eliot Courtney.
#ifndef COMMON_TEST_H_
#define COMMON_TEST_H_

#include <cstdint>

#include "compute/energy/model.h"

namespace mrna {

inline constexpr uint32_t T04_P1_MODEL_HASH = 0x768ab8e1;
inline constexpr uint32_t T04_P2_MODEL_HASH = 0x3bc125da;

// Make sure to use Range(0, NUM_TEST_MODELS) if making a parameterised test
// with all models in g_em, since g_em is initialized at runtime.
inline constexpr int NUM_TEST_MODELS = 5;
inline energy::EnergyModelPtr g_em[NUM_TEST_MODELS];

#if ENERGY_PRECISION == 1
inline energy::EnergyModelPtr t04_p1;
#elif ENERGY_PRECISION == 2
inline energy::EnergyModelPtr t04_p2;
#endif

}  // namespace mrna

#endif  // COMMON_TEST_H_
