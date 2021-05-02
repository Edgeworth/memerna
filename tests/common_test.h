// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#ifndef MEMERNA_COMMON_TEST_H
#define MEMERNA_COMMON_TEST_H

#include <cstdint>
#include <iostream>

#include "energy/energy_model.h"
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

namespace memerna {

const uint32_t T04_MODEL_HASH = 0xcf2538f8;
extern energy::EnergyModelPtr g_em;
extern std::vector<energy::EnergyModelPtr> g_ems;

std::ostream& operator<<(std::ostream& os, const secondary_t& s);
std::ostream& operator<<(std::ostream& os, const computed_t& computed);
}  // namespace memerna

#endif  // MEMERNA_COMMON_TEST_H
