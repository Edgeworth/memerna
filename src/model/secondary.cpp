// Copyright 2022 E.
#include "model/secondary.h"

#include <stack>

#include "model/primary.h"
#include "util/error.h"

namespace mrna {

Secondary Secondary::FromDotBracket(const std::string& pairs_str) {
  Secondary s(pairs_str.size());
  std::stack<int> stk;
  for (int i = 0; i < static_cast<int>(pairs_str.size()); ++i) {
    if (pairs_str[i] == '(') {
      stk.push(i);
    } else if (pairs_str[i] == ')') {
      verify(!stk.empty(), "unmatched bracket");
      s[i] = stk.top();
      s[stk.top()] = i;
      stk.pop();
    }
  }
  return s;
}

std::string Secondary::ToDotBracket() const {
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

std::tuple<Primary, Secondary> ParsePrimaryDotBracket(
    const std::string& prim_str, const std::string& pairs_str) {
  verify(prim_str.size() == pairs_str.size(), "requires rna length to be the same as pairs length");
  return {Primary::FromString(prim_str), Secondary::FromDotBracket(pairs_str)};
}

}  // namespace mrna
