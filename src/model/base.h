// Copyright 2016 Eliot Courtney.
#ifndef MODEL_BASE_H_
#define MODEL_BASE_H_

#include <cstdint>

namespace mrna {

using Base = int8_t;

// Don't ever change these values.
const Base A = 0, C = 1, G = 2, U = 3;
const int A_b = 1 << A, C_b = 1 << C, G_b = 1 << G, U_b = 1 << U;

inline bool IsBase(Base b) { return b >= 0 && b <= 3; }

inline bool IsPairOf(Base a, Base b, int basesA, int basesB) {
  return (basesA & (1 << a)) && (basesB & 1 << b);
}

inline bool IsUnorderedOf(Base a, Base b, int basesA, int basesB) {
  return IsPairOf(a, b, basesA, basesB) || IsPairOf(b, a, basesA, basesB);
}

inline bool CanPair(Base a, Base b) {
  int combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (G_b | C_b) || combined == (A_b | U_b);
}

inline bool IsAuGu(Base a, Base b) {
  int combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (A_b | U_b);
}

inline bool IsWatsonCrick(Base a, Base b) {
  int combined = (1 << a) | (1 << b);
  return combined == (A_b | U_b) || combined == (G_b | C_b);
}

inline bool IsGu(Base a, Base b) { return (a == G && b == U) || (a == U && b == G); }

Base CharToBase(char c);

char BaseToChar(Base b);

}  // namespace mrna

#endif  // MODEL_BASE_H_
