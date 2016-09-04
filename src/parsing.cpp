#include "parsing.h"

namespace memerna {
namespace parsing {

using namespace energy;

primary_t StringToPrimary(const std::string& s) {
  primary_t r(s.size());
  for (int i = 0; i < int(s.size()); ++i) {
    r[i] = CharToBase(s[i]);
    verify_expr(r[i] != -1, "unexpected base %c", s[i]);
  }
  return r;
}

secondary_t ParseDotBracketSecondary(const std::string& prim_str, const std::string& pairs_str) {
  verify_expr(prim_str.size() == pairs_str.size(), "requires rna length to be the same as pairs length");
  return {StringToPrimary(prim_str), DotBracketToPairs(pairs_str)};
}

std::vector<int> DotBracketToPairs(const std::string& pairs_str) {
  std::vector<int> pairs(pairs_str.size(), -1);
  std::stack<int> s;
  for (int i = 0; i < int(pairs_str.size()); ++i) {
    if (pairs_str[i] == '(') {
      s.push(i);
    } else if (pairs_str[i] == ')') {
      verify_expr(!s.empty(), "unmatched bracket");
      pairs[i] = s.top();
      pairs[s.top()] = i;
      s.pop();
    }
  }
  return pairs;
}

std::string PairsToDotBracket(const std::vector<int>& pairs) {
  std::string s(pairs.size(), '.');
  for (int i = 0; i < int(pairs.size()); ++i) {
    if (pairs[i] == -1) continue;
    if (pairs[i] < i) s[i] = ')';
    else s[i] = '(';
  }
  return s;
}

std::string ComputedToCtdString(const computed_t& computed) {
  const auto& p = computed.s.p;
  std::string s = PairsToDotBracket(p);
  for (int i = 0; i < int(p.size()); ++i) {
    switch (computed.base_ctds[i]) {
      case CTD_NA:
      case CTD_UNUSED:
        break;
      case CTD_3_DANGLE:
        s[p[i] + 1] = '3';
        break;
      case CTD_5_DANGLE:
        s[i - 1] = '5';
        break;
      case CTD_TERMINAL_MISMATCH:
        s[i - 1] = 'm';
        s[p[i] + 1] = 'm';
        break;
      case CTD_LEFT_MISMATCH_COAX_WITH_NEXT:
        s[i] = '>';
        s[i - 1] = 'm';
        s[p[i] - 1] = 'm';
        break;
      case CTD_LEFT_MISMATCH_COAX_WITH_PREV:
        s[i] = '<';
        break;
      case CTD_RIGHT_MISMATCH_COAX_WITH_NEXT:
        s[i] = '>';
        break;
      case CTD_RIGHT_MISMATCH_COAX_WITH_PREV:
        s[i] = '<';
        s[i - 1] = 'm';
        s[p[i] - 1] = 'm';
        break;
      case CTD_FLUSH_COAX_WITH_NEXT:
        s[i] = '>';
        break;
      case CTD_FLUSH_COAX_WITH_PREV:
        s[i] = '<';
        break;
      default:
        verify_expr(false, "bug");
    }
  }
  return s;
}

std::string PrimaryToString(const primary_t& r) {
  std::string s;
  s.resize(r.size());
  for (int i = 0; i < int(r.size()); ++i) {
    s[i] = BaseToChar(r[i]);
  }
  return s;
}


}
}
