// Copyright 2022 E.
#include "model/primary.h"

#include <random>

#include "util/error.h"

namespace mrna {

Primary Primary::Random(int length) {
  static thread_local std::mt19937 eng;
  std::uniform_int_distribution<int> dist(0, 3);
  Primary r(length);
  for (int i = 0; i < length; ++i) r[i] = Base(dist(eng));
  return r;
}

Primary Primary::FromString(const std::string& s) {
  Primary r(s.size());
  for (int i = 0; i < static_cast<int>(s.size()); ++i) {
    r[i] = CharToBase(s[i]);
    verify(r[i] != -1, "unexpected base %c", s[i]);
  }
  return r;
}

std::string Primary::ToString() const {
  std::string s;
  s.resize(size());
  for (int i = 0; i < static_cast<int>(size()); ++i) s[i] = BaseToChar(data_[i]);
  return s;
}

}  // namespace mrna
