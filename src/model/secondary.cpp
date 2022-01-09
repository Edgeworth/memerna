// Copyright 2022 Eliot Courtney.
#include "model/secondary.h"

#include <stack>

#include "util/macros.h"

namespace mrna {

Secondary ParseDotBracketSecondary(const std::string& prim_str, const std::string& pairs_str) {
  verify(prim_str.size() == pairs_str.size(), "requires rna length to be the same as pairs length");
  return {StringToPrimary(prim_str), DotBracketToPairs(pairs_str)};
}

std::vector<int> DotBracketToPairs(const std::string& pairs_str) {
  std::vector<int> pairs(pairs_str.size(), -1);
  std::stack<int> s;
  for (int i = 0; i < static_cast<int>(pairs_str.size()); ++i) {
    if (pairs_str[i] == '(') {
      s.push(i);
    } else if (pairs_str[i] == ')') {
      verify(!s.empty(), "unmatched bracket");
      pairs[i] = s.top();
      pairs[s.top()] = i;
      s.pop();
    }
  }
  return pairs;
}

std::string PairsToDotBracket(const std::vector<int>& pairs) {
  std::string s(pairs.size(), '.');
  for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
    if (pairs[i] == -1) continue;
    if (pairs[i] < i)
      s[i] = ')';
    else
      s[i] = '(';
  }
  return s;
}

}  // namespace mrna
