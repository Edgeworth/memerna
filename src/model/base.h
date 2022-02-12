// Copyright 2016 Eliot Courtney.
#ifndef MODEL_BASE_H_
#define MODEL_BASE_H_

#include <cassert>
#include <cstdint>
#include <limits>

#include "util/error.h"

namespace mrna {

using Base = uint8_t;
using BaseMask = uint32_t;

// Don't ever change these values.
constexpr Base A = 0, C = 1, G = 2, U = 3, INVALID_BASE = std::numeric_limits<uint8_t>::max();

constexpr BaseMask A_b = 1U << A, C_b = 1U << C, G_b = 1U << G, U_b = 1U << U;

inline constexpr bool IsBase(Base b) { return b >= 0 && b <= 3; }

inline constexpr bool IsPairOf(Base a, Base b, BaseMask basesA, BaseMask basesB) {
  return (basesA & (1U << a)) && (basesB & (1U << b));
}

inline constexpr bool IsUnorderedOf(Base a, Base b, BaseMask basesA, BaseMask basesB) {
  return IsPairOf(a, b, basesA, basesB) || IsPairOf(b, a, basesA, basesB);
}

inline constexpr bool IsPair(Base a, Base b) {
  int combined = (1U << a) | (1U << b);
  return combined == (G_b | U_b) || combined == (G_b | C_b) || combined == (A_b | U_b);
}

inline constexpr bool IsAuGuPair(Base a, Base b) {
  int combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (A_b | U_b);
}

inline constexpr bool IsGuPair(Base a, Base b) { return (a == G && b == U) || (a == U && b == G); }

inline constexpr bool IsWcPair(Base a, Base b) {
  int combined = (1U << a) | (1U << b);
  return combined == (A_b | U_b) || combined == (G_b | C_b);
}

inline constexpr bool IsGu(Base b) { return b == G || b == U; }

// Returns the Watson-Crick complement of the given base.
inline constexpr Base WcPair(Base b) { return b ^ 3U; }

// Returns the GU complement of the given base.
inline constexpr Base GuPair(Base b) {
  assert(IsGu(b));
  return b ^ 1U;
}

Base CharToBase(char c);

char BaseToChar(Base b);

}  // namespace mrna

#endif  // MODEL_BASE_H_
