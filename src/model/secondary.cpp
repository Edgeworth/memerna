// Copyright 2022 Eliot Courtney.
#include "model/secondary.h"

#include <string>
#include <tuple>
#include <vector>

#include "model/base.h"
#include "model/constants.h"
#include "model/primary.h"
#include "util/error.h"
#include "util/util.h"

namespace mrna {

void Secondary::Verify(const Primary& r) const {
  verify(r.size() == size(), "sequence and secondary structure must be the same length");
  for (int i = 0; i < size(); ++i) {
    if (!IsOpeningPair(i)) continue;
    const int j = data_[i];
    verify(IsPair(r[i], r[j]), "non-canonical pair ({}, {})", i, j);
    verify(j - i - 1 >= HAIRPIN_MIN_SZ, "pair ({}, {}) encloses fewer than {} bases", i, j,
        HAIRPIN_MIN_SZ);
  }
}

Secondary Secondary::FromDb(const std::string& pairs_str) {
  Secondary s(pairs_str.size());
  const int N = s.size();
  std::vector<int> stk;
  for (int i = 0; i < N; ++i) {
    if (pairs_str[i] == '(') {
      stk.push_back(i);
    } else if (pairs_str[i] == ')') {
      verify(!stk.empty(), "unmatched closing bracket at position {}", i);
      s[i] = As<Index>(stk.back());
      s[stk.back()] = As<Index>(i);
      stk.pop_back();
    } else {
      verify(pairs_str[i] == '.', "unexpected character '{}' at position {}", pairs_str[i], i);
    }
  }
  verify(stk.empty(), "unmatched opening bracket at position {}", stk.back());
  return s;
}

std::string Secondary::ToDb() const {
  std::string db(size(), '.');
  for (int i = 0; i < size(); ++i) {
    if (IsUnpaired(i)) continue;
    db[i] = IsClosingPair(i) ? ')' : '(';
  }
  return db;
}

std::tuple<Primary, Secondary> ParseSeqDb(
    const std::string& prim_str, const std::string& pairs_str) {
  verify(prim_str.size() == pairs_str.size(), "requires rna length to be the same as pairs length");
  return {Primary::FromSeq(prim_str), Secondary::FromDb(pairs_str)};
}

}  // namespace mrna
