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
#ifndef MEMERNA_BASE_H
#define MEMERNA_BASE_H

#include <cstdint>
#include <cstring>
#include <memory>
#include <random>

#include "common.h"

namespace memerna {

// Don't ever change these values.
const base_t A = 0, C = 1, G = 2, U = 3;
const int A_b = 1 << A, C_b = 1 << C, G_b = 1 << G, U_b = 1 << U;

inline bool IsBase(base_t b) {
  return b >= 0 && b <= 3;
}

inline bool IsPairOf(base_t a, base_t b, int basesA, int basesB) {
  return (basesA & (1 << a)) && (basesB & 1 << b);
}

inline bool IsUnorderedOf(base_t a, base_t b, int basesA, int basesB) {
  return IsPairOf(a, b, basesA, basesB) || IsPairOf(b, a, basesA, basesB);
}

inline bool CanPair(base_t a, base_t b) {
  int combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (G_b | C_b) || combined == (A_b | U_b);
}

inline bool IsAuGu(base_t a, base_t b) {
  int combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (A_b | U_b);
}

inline bool IsWatsonCrick(base_t a, base_t b) {
  int combined = (1 << a) | (1 << b);
  return combined == (A_b | U_b) || combined == (G_b | C_b);
}

inline bool IsGu(base_t a, base_t b) { return (a == G && b == U) || (a == U && b == G); }

base_t CharToBase(char c);

char BaseToChar(base_t b);

template <typename RandomEngine>
primary_t GenerateRandomPrimary(int length, RandomEngine& eng) {
  std::uniform_int_distribution<int> dist(0, 3);
  primary_t r(std::size_t(length), 0);
  for (int i = 0; i < length; ++i)
    r[i] = base_t(dist(eng));
  return r;
}
}

#endif  // MEMERNA_BASE_H
