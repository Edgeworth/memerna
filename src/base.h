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

template <typename RandomEngine> primary_t GenerateRandomPrimary(int length, RandomEngine& eng) {
  std::uniform_int_distribution<int> dist(0, 3);
  primary_t r(std::size_t(length), 0);
  for (int i = 0; i < length; ++i) r[i] = base_t(dist(eng));
  return r;
}
}

#endif  // MEMERNA_BASE_H
