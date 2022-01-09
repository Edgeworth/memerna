// Copyright 2022 Eliot Courtney.
#include "model/secondary.h"

#include <stack>

#include "util/macros.h"

namespace mrna {

std::tuple<Primary, Secondary> ParsePrimaryDotBracket(
    const std::string& prim_str, const std::string& pairs_str) {
  verify(prim_str.size() == pairs_str.size(), "requires rna length to be the same as pairs length");
  return {StringToPrimary(prim_str), DotBracketToSecondary(pairs_str)};
}

Secondary DotBracketToSecondary(const std::string& pairs_str) {
  Secondary s(pairs_str.size(), -1);
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

std::string SecondaryToDotBracket(const Secondary& s) {
  std::string db(s.size(), '.');
  for (int i = 0; i < static_cast<int>(s.size()); ++i) {
    if (s[i] == -1) continue;
    if (s[i] < i)
      db[i] = ')';
    else
      db[i] = '(';
  }
  return db;
}

}  // namespace mrna
