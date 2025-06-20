// Copyright 2022 E.
#include "model/primary.h"

#include <optional>
#include <random>
#include <string>

#include "util/error.h"

namespace mrna {

Primary Primary::Random(int length, std::mt19937& eng) {
  std::uniform_int_distribution<int> dist(0, 3);
  Primary r(length);
  for (int i = 0; i < length; ++i) r[i] = Base(dist(eng));
  return r;
}

Primary Primary::FromSeq(const std::string& s) {
  Primary r(s.size());
  for (int i = 0; i < static_cast<int>(s.size()); ++i) {
    const auto base = CharToBase(s[i]);
    verify(base.has_value(), "unexpected base {}", s[i]);
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

void Primary::Increment() {
  bool carry = true;
  for (auto& base : data_) {
    base++;
    if (base != MAX_BASE) {
      carry = false;
      break;
    }
    base = MIN_BASE;
  }
  if (carry) data_.push_back(0);
}

}  // namespace mrna
