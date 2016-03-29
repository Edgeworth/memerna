#include "base.h"
#include "parsing.h"

namespace memerna {

void Init() {
  parsing::ParseStackingEnergiesFromFile("data/stacking.data");
}

base_t BaseFromChar(char c) {
  switch (c) {
    case 'A':
      return A;
    case 'C':
      return C;
    case 'G':
      return G;
    case 'U':
      return U;
    default:
      return -1;
  }
}

char BaseToChar(base_t b) {
  if (b < 0 || b > 4) return '?';
  return "ACGU"[b];
}

}
