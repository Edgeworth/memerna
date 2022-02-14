// Copyright 2022 Eliot Courtney.
#include "model/primary.h"

#include <optional>
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

Primary Primary::FromSeq(const std::string& s) {
  Primary r(s.size());
  for (int i = 0; i < static_cast<int>(s.size()); ++i) {
    const auto base = CharToBase(s[i]);
    verify(base.has_value(), "unexpected base %c", s[i]);
    r[i] = *base;
  }
  return r;
}

std::string Primary::ToSeq() const {
  std::string s;
  s.resize(size());
  for (int i = 0; i < static_cast<int>(size()); ++i) s[i] = BaseToChar(data_[i]).value();
  return s;
}

}  // namespace mrna
