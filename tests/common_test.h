// Copyright 2016 Eliot Courtney.
#ifndef TESTS_COMMON_TEST_H_
#define TESTS_COMMON_TEST_H_

#include <cstdint>
#include <iostream>
#include <vector>

#include "compute/energy/model.h"
#include "gtest/gtest.h"

#define ONLY_FOR_THIS_MODEL(em_, hash_)                                                          \
  do {                                                                                           \
    auto our_hash = em_->Checksum();                                                             \
    if (our_hash != hash_) {                                                                     \
      printf("Skipping energy model specific tests: %#010x != " #hash_ " (%#010x).\n", our_hash, \
          hash_);                                                                                \
      return;                                                                                    \
    }                                                                                            \
  } while (0)

namespace mrna {

const uint32_t T04_MODEL_HASH = 0x9eabeccf;
extern energy::EnergyModelPtr g_em;

std::ostream& operator<<(std::ostream& os, const Secondary& s);
std::ostream& operator<<(std::ostream& os, const Computed& computed);

}  // namespace mrna

#endif  // TESTS_COMMON_TEST_H_
