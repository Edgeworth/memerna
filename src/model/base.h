// Copyright 2016 Eliot Courtney.
#ifndef MODEL_BASE_H_
#define MODEL_BASE_H_

#include <cstdint>
#include <cstring>
#include <memory>
#include <random>

namespace mrna {

typedef int8_t base_t;

// Don't ever change these values.
const base_t A = 0, C = 1, G = 2, U = 3;
const int A_b = 1 << A, C_b = 1 << C, G_b = 1 << G, U_b = 1 << U;

inline bool IsBase(base_t b) { return b >= 0 && b <= 3; }

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

}  // namespace mrna

#endif  // MODEL_BASE_H_
