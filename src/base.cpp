// Copyright 2016 E.
#include "base.h"

namespace mrna {

base_t CharToBase(char c) {
  switch (c) {
  case 'A': return A;
  case 'C': return C;
  case 'G': return G;
  case 'U': return U;
  default: return -1;
  }
}

char BaseToChar(base_t b) {
  if (b < 0 || b > 4) return '?';
  return "ACGU"[b];
}

}  // namespace mrna
