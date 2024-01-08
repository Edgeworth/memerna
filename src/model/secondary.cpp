// Copyright 2022 Eliot Courtney.
#include "model/secondary.h"

#include <vector>

#include "model/primary.h"
#include "util/error.h"

namespace mrna {

Secondary Secondary::FromDb(const std::string& pairs_str) {
  Secondary s(pairs_str.size());
  std::vector<int> stk;
  for (int i = 0; i < static_cast<int>(pairs_str.size()); ++i) {
    if (pairs_str[i] == '(') {
      stk.push_back(i);
    } else if (pairs_str[i] == ')') {
      verify(!stk.empty(), "unmatched bracket");
      s[i] = stk.back();
      s[stk.back()] = i;
      stk.pop_back();
    }
  }
  return s;
}

std::string Secondary::ToDb() const {
  std::string db(size(), '.');
  for (int i = 0; i < static_cast<int>(size()); ++i) {
    if (data_[i] == -1) continue;
    if (data_[i] < i)
      db[i] = ')';
    else
      db[i] = '(';
  }
  return db;
}

std::tuple<Primary, Secondary> ParseSeqDb(
    const std::string& prim_str, const std::string& pairs_str) {
  verify(prim_str.size() == pairs_str.size(), "requires rna length to be the same as pairs length");
  return {Primary::FromSeq(prim_str), Secondary::FromDb(pairs_str)};
}

}  // namespace mrna
