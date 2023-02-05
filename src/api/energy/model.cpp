// Copyright 2022 Eliot Courtney.
#include "api/energy/model.h"

#include <ostream>
#include <string>

#include "util/string.h"

namespace mrna::erg {

std::istream& operator>>(std::istream& is, ModelKind& kind) {
  std::string s;
  is >> s;
  s = ToLower(s);
  if (s == "t04p1" || s == "t04p2" || s == "t12p2")
    kind = ModelKind::T04_LIKE;
  else if (s == "t22p2")
    kind = ModelKind::T22_LIKE;
  else
    error("Invalid energy model {}", s);
  return is;
}

std::ostream& operator<<(std::ostream& os, ModelKind kind) {
  switch (kind) {
  case ModelKind::T04_LIKE: return os << "t04-like";
  case ModelKind::T22_LIKE: return os << "t22-like";
  default: bug();
  }
}

}  // namespace mrna::erg
