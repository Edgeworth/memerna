#ifndef MEMERNA_BASE_H
#define MEMERNA_BASE_H

#include <cstdint>
#include <memory>
#include <cstring>

#include "common.h"

namespace memerna {

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

base_t CharToBase(char c);

char BaseToChar(base_t b);

}

#endif //MEMERNA_BASE_H
