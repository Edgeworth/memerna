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
  verify_expr(
      prim_str.size() == pairs_str.size(), "requires rna length to be the same as pairs length");
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
    if (pairs[i] < i)
      s[i] = ')';
    else
      s[i] = '(';
  }
  return s;
}

std::string ComputedToCtdString(const computed_t& computed) {
  const auto& p = computed.s.p;
  std::string s = PairsToDotBracket(p);
  for (int i = 0; i < int(p.size()); ++i) {
    bool closing = false;
    if (s[i] == '(') s[i] = '[';
    if (s[i] == ')') {
      s[i] = ']';
      closing = true;
    }
    switch (computed.base_ctds[i]) {
    case CTD_NA:
    case CTD_UNUSED: break;
    case CTD_3_DANGLE: s[p[i] + 1] = '3'; break;
    case CTD_5_DANGLE: s[i - 1] = '5'; break;
    case CTD_MISMATCH:
      s[i - 1] = 'm';
      s[p[i] + 1] = 'M';
      break;
    case CTD_LCOAX_WITH_NEXT:
      s[i] = closing ? 'N' : 'n';
      s[i - 1] = 'm';
      s[p[i] + 1] = 'M';
      break;
    case CTD_LCOAX_WITH_PREV: s[i] = closing ? 'P' : 'p'; break;
    case CTD_RCOAX_WITH_NEXT: s[i] = closing ? 'N' : 'n'; break;
    case CTD_RCOAX_WITH_PREV:
      s[i] = closing ? 'P' : 'p';
      s[i - 1] = 'm';
      s[p[i] + 1] = 'M';
      break;
    case CTD_FCOAX_WITH_NEXT: s[i] = closing ? 'N' : 'n'; break;
    case CTD_FCOAX_WITH_PREV: s[i] = closing ? 'P' : 'p'; break;
    default: verify_expr(false, "bug");
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

// Breaking ambiguous case for the old format:
// >(   )>   )<<   )(   )>   )<   )
// (m(   )m(   )m(   )m)

computed_t ParseCtdComputed(const std::string& prim_str, const std::string& pairs_str) {
  computed_t computed(StringToPrimary(prim_str));
  verify_expr(
      prim_str.size() == pairs_str.size(), "primary and pairs string need to be the same length");
  const int N = int(prim_str.size());
  std::stack<int> s;
  auto& p = computed.s.p;
  for (int i = 0; i < N; ++i) {
    char c = pairs_str[i];
    if (c == 'p' || c == 'n' || c == '[') {
      s.push(i);
    } else if (c == 'P' || c == 'N' || c == ']') {
      verify_expr(!s.empty(), "invalid input");
      p[s.top()] = i;
      p[i] = s.top();
      s.pop();
    }
  }
  verify_expr(s.empty(), "invalid input");

  const std::string allowed_characters = "[]nNpP.mM35";
  const std::string branch_characters = "nNpP[]";
  for (int i = 0; i < N; ++i) {
    char c = pairs_str[i];
    // Coaxial stacking handling.
    if (c == 'P' || c == 'p') {
      int prev = i - 1;
      bool mismatch = false;
      verify_expr(prev >= 0, "invalid input");
      if (p[prev] == -1) {
        --prev;
        mismatch = true;
      }
      verify_expr(
          prev >= 0 && p[prev] != -1 && (pairs_str[p[prev]] == 'n' || pairs_str[p[prev]] == 'N'),
          "invalid input");
      verify_expr(computed.base_ctds[i] == CTD_NA && computed.base_ctds[p[prev]] == CTD_NA,
          "invalid input");
      bool left_mismatch = p[prev] - 1 >= 0 && pairs_str[p[prev] - 1] == 'm';
      bool right_mismatch = p[i] + 1 < N && pairs_str[p[i] + 1] == 'M';
      verify_expr(!mismatch || (left_mismatch ^ right_mismatch), "invalid input");
      if (mismatch) {
        if (left_mismatch) {
          computed.base_ctds[i] = CTD_LCOAX_WITH_PREV;
          computed.base_ctds[p[prev]] = CTD_LCOAX_WITH_NEXT;
        } else {
          computed.base_ctds[i] = CTD_RCOAX_WITH_PREV;
          computed.base_ctds[p[prev]] = CTD_RCOAX_WITH_NEXT;
        }
      } else {
        computed.base_ctds[i] = CTD_FCOAX_WITH_PREV;
        computed.base_ctds[p[prev]] = CTD_FCOAX_WITH_NEXT;
      }
    }
    if (c == '3') {
      // If we optimistically set an exterior loop to CTD_UNUSED, we might want to rewrite it here.
      verify_expr(i - 1 >= 0 && p[i - 1] != -1 && (computed.base_ctds[p[i - 1]] == CTD_NA ||
                                                      computed.base_ctds[p[i - 1]] == CTD_UNUSED),
          "invalid input");
      computed.base_ctds[p[i - 1]] = CTD_3_DANGLE;
    }
    if (c == '5') {
      verify_expr(i + 1 < N && computed.base_ctds[i + 1] == CTD_NA, "invalid input");
      computed.base_ctds[i + 1] = CTD_5_DANGLE;
    }
    // Opening a branch.
    if (c == 'p' || c == 'n' || c == '[') {
      if (!s.empty())
        ++s.top();  // One more branch here.
      else if (computed.base_ctds[i] == CTD_NA && c == '[')
        computed.base_ctds[i] = CTD_UNUSED;
      s.push(0);  // Number of branches for this part.
    }
    // Closing a branch.
    if (c == 'P' || c == 'N' || c == ']') {
      // Set explicitly unused CTDs for any child branches if this is a multiloop.
      // Also this branch as an outer loop.
      if (s.top() >= 2) {
        if (computed.base_ctds[i] == CTD_NA) computed.base_ctds[i] = CTD_UNUSED;
        for (int j = p[i] + 1; j < i; ++j) {
          if (p[j] == -1) continue;
          if (computed.base_ctds[j] == CTD_NA) computed.base_ctds[j] = CTD_UNUSED;
          j = p[j];
        }
      }
      s.pop();
    }
    if (allowed_characters.find(c) == std::string::npos)
      verify_expr(false, "invalid input '%c'", pairs_str[i]);
  }
  // Add in the terminal mismatches.
  for (int i = 0; i < N; ++i) {
    char c = pairs_str[i];
    if (c == 'm') {
      bool right_branch_viable =
          i + 1 < N && p[i + 1] != -1 && p[i + 1] + 1 < N && pairs_str[p[i + 1] + 1] == 'M';
      verify_expr(right_branch_viable && computed.base_ctds[i + 1] != CTD_NA, "invalid input");
      if (computed.base_ctds[i + 1] == CTD_UNUSED) computed.base_ctds[i + 1] = CTD_MISMATCH;
    }
  }
  return computed;
}

bool IsCtdString(const std::string& pairs_str) {
  for (auto c : pairs_str)
    if (c == '(') return false;
  return true;
}
}
}
