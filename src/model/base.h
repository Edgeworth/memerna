// Copyright 2016 E.
#ifndef MODEL_BASE_H_
#define MODEL_BASE_H_

#include <cassert>
#include <cstdint>

#include "util/error.h"

namespace mrna {

using Base = int8_t;

// Don't ever change these values.
constexpr Base A = 0, C = 1, G = 2, U = 3;
constexpr int A_b = 1 << A, C_b = 1 << C, G_b = 1 << G, U_b = 1 << U;

inline constexpr bool IsBase(Base b) { return b >= 0 && b <= 3; }

inline constexpr bool IsPairOf(Base a, Base b, int basesA, int basesB) {
  return (basesA & (1 << a)) && (basesB & (1 << b));
}

inline constexpr bool IsUnorderedOf(Base a, Base b, int basesA, int basesB) {
  return IsPairOf(a, b, basesA, basesB) || IsPairOf(b, a, basesA, basesB);
}

inline constexpr bool IsPair(Base a, Base b) {
  int combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (G_b | C_b) || combined == (A_b | U_b);
}

inline constexpr bool IsAuGuPair(Base a, Base b) {
  int combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (A_b | U_b);
}

inline constexpr bool IsGuPair(Base a, Base b) { return (a == G && b == U) || (a == U && b == G); }

inline constexpr bool IsWcPair(Base a, Base b) {
  int combined = (1 << a) | (1 << b);
  return combined == (A_b | U_b) || combined == (G_b | C_b);
}

inline constexpr bool IsGu(Base b) { return b == G || b == U; }

// Returns the Watson-Crick complement of the given base.
inline constexpr Base WcPair(Base b) { return b ^ 3; }

// Returns the GU complement of the given base.
inline constexpr Base GuPair(Base b) {
  assert(IsGu(b));
  return b ^ 1;
}

Base CharToBase(char c);

char BaseToChar(Base b);

}  // namespace mrna

#endif  // MODEL_BASE_H_
