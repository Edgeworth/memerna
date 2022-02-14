// Copyright 2016 Eliot Courtney.
#include "model/base.h"

namespace mrna {

std::optional<Base> CharToBase(char c) {
  switch (c) {
  case 'A': return A;
  case 'C': return C;
  case 'G': return G;
  case 'U': return U;
  default: return std::nullopt;
  }
}

std::optional<char> BaseToChar(Base b) {
  switch (b) {
  case A: return 'A';
  case C: return 'C';
  case G: return 'G';
  case U: return 'U';
  default: return std::nullopt;
  }
}

}  // namespace mrna
