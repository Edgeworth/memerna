// Copyright 2016 Eliot Courtney.
#ifndef TESTS_COMMON_TEST_H_
#define TESTS_COMMON_TEST_H_

#include <cstdint>

#include "compute/energy/model.h"

#define ONLY_FOR_THIS_MODEL(em_, hash_)                                                          \
  do {                                                                                           \
    auto our_hash = (em_).Checksum();                                                            \
    if (our_hash != (hash_)) {                                                                   \
      printf("Skipping energy model specific tests: %#010x != " #hash_ " (%#010x).\n", our_hash, \
          hash_);                                                                                \
      return;                                                                                    \
    }                                                                                            \
  } while (0)

namespace mrna {

inline constexpr uint32_t T04_MODEL_HASH = 0x9eabeccf;

// Make sure to use Range(0, NUM_TEST_MODELS) if making a parameterised test
// with all models in g_em, since g_em is initialized at runtime.
inline constexpr int NUM_TEST_MODELS = 5;
inline energy::EnergyModelPtr g_em[NUM_TEST_MODELS];
inline energy::EnergyModelPtr t04;

}  // namespace mrna

#endif  // TESTS_COMMON_TEST_H_
