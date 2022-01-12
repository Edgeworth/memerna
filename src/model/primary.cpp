// Copyright 2022 Eliot Courtney.
#include "model/primary.h"

#include "util/macros.h"

namespace mrna {

Primary GenerateRandomPrimary(int length) {
  static thread_local std::mt19937 eng;
  std::uniform_int_distribution<int> dist(0, 3);
  Primary r((std::size_t(length)));
  for (int i = 0; i < length; ++i) r[i] = Base(dist(eng));
  return r;
}

Primary StringToPrimary(const std::string& s) {
  Primary r(s.size());
  for (int i = 0; i < static_cast<int>(s.size()); ++i) {
    r[i] = CharToBase(s[i]);
    verify(r[i] != -1, "unexpected base %c", s[i]);
  }
  return r;
}

std::string PrimaryToString(const Primary& r) {
  std::string s;
  s.resize(r.size());
  for (int i = 0; i < static_cast<int>(r.size()); ++i) { s[i] = BaseToChar(r[i]); }
  return s;
}

}  // namespace mrna
