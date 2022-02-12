// Copyright 2016 Eliot Courtney.
#include "model/base.h"

#include <limits>

namespace mrna {

Base CharToBase(char c) {
  switch (c) {
  case 'A': return A;
  case 'C': return C;
  case 'G': return G;
  case 'U': return U;
  default: return INVALID_BASE;
  }
}

char BaseToChar(Base b) {
  if (b < 0 || b > 4) return '?';
  return "ACGU"[b];
}

}  // namespace mrna
