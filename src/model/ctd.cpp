// Copyright 2022 Eliot Courtney.
#include "model/ctd.h"

#include <functional>
#include <stack>
#include <utility>

#include "model/primary.h"
#include "util/error.h"

namespace mrna {

std::string Ctds::ToString(const Secondary& s) const {
  std::string str(s.size(), '.');
  for (int i = 0; i < static_cast<int>(s.size()); ++i) {
    if (s[i] == -1) continue;
    const bool closing = s[i] < i;
    if (closing)
      str[i] = ']';
    else
      str[i] = '[';
    switch (data_[i]) {
    case CTD_NA:
    case CTD_UNUSED: break;
    case CTD_3_DANGLE: str[s[i] + 1] = '3'; break;
    case CTD_5_DANGLE: str[i - 1] = '5'; break;
    case CTD_MISMATCH:
      str[i - 1] = 'm';
      str[s[i] + 1] = 'M';
      break;
    case CTD_LCOAX_WITH_NEXT:
      str[i] = closing ? 'N' : 'n';
      str[i - 1] = 'm';
      str[s[i] + 1] = 'M';
      break;
    case CTD_LCOAX_WITH_PREV: str[i] = closing ? 'P' : 'p'; break;
    case CTD_RC_WITH_NEXT: str[i] = closing ? 'N' : 'n'; break;
    case CTD_RC_WITH_PREV:
      str[i] = closing ? 'P' : 'p';
      str[i - 1] = 'm';
      str[s[i] + 1] = 'M';
      break;
    case CTD_FCOAX_WITH_NEXT: str[i] = closing ? 'N' : 'n'; break;
    case CTD_FCOAX_WITH_PREV: str[i] = closing ? 'P' : 'p'; break;
    default: bug();
    }
  }
  return str;
}

bool Ctds::IsCtdString(const std::string& ctd_str) {
  return std::ranges::none_of(ctd_str, [](char c) { return c == '('; });
}

// Breaking ambiguous case for the old format:
// >(   )>   )<<   )(   )>   )<   )
// (m(   )m(   )m(   )m)

std::tuple<Primary, Secondary, Ctds> ParseSeqCtdString(
    const std::string& prim_str, const std::string& ctd_str) {
  auto r = Primary::FromSeq(prim_str);
  verify(prim_str.size() == ctd_str.size(), "primary and pairs string need to be the same length");
  const int N = static_cast<int>(r.size());
  std::stack<int> stk;
  Secondary s(N);
  for (int i = 0; i < N; ++i) {
    const char c = ctd_str[i];
    if (c == 'p' || c == 'n' || c == '[') {
      stk.push(i);
    } else if (c == 'P' || c == 'N' || c == ']') {
      verify(!stk.empty(), "invalid input");
      s[stk.top()] = i;
      s[i] = stk.top();
      stk.pop();
    }
  }
  verify(stk.empty(), "invalid input");

  const std::string allowed_characters = "[]nNpP.mM35";
  const std::string branch_characters = "nNpP[]";
  Ctds ctd(N);
  for (int i = 0; i < N; ++i) {
    const char c = ctd_str[i];
    // Coaxial stacking handling.
    if (c == 'P' || c == 'p') {
      int prev = i - 1;
      bool mismatch = false;
      verify(prev >= 0, "invalid input");
      if (s[prev] == -1) {
        --prev;
        mismatch = true;
      }
      verify(prev >= 0 && s[prev] != -1 && (ctd_str[s[prev]] == 'n' || ctd_str[s[prev]] == 'N'),
          "invalid input");
      verify(ctd[i] == CTD_NA && ctd[s[prev]] == CTD_NA, "invalid input");
      const bool left_mismatch = s[prev] - 1 >= 0 && ctd_str[s[prev] - 1] == 'm';
      const bool right_mismatch = s[i] + 1 < N && ctd_str[s[i] + 1] == 'M';
      verify(!mismatch || (left_mismatch ^ right_mismatch), "invalid input");
      if (mismatch) {
        if (left_mismatch) {
          ctd[i] = CTD_LCOAX_WITH_PREV;
          ctd[s[prev]] = CTD_LCOAX_WITH_NEXT;
        } else {
          ctd[i] = CTD_RC_WITH_PREV;
          ctd[s[prev]] = CTD_RC_WITH_NEXT;
        }
      } else {
        ctd[i] = CTD_FCOAX_WITH_PREV;
        ctd[s[prev]] = CTD_FCOAX_WITH_NEXT;
      }
    }
    if (c == '3') {
      // If we optimistically set an exterior loop to CTD_UNUSED, we might want to rewrite it here.
      verify(
          i - 1 >= 0 && s[i - 1] != -1 && (ctd[s[i - 1]] == CTD_NA || ctd[s[i - 1]] == CTD_UNUSED),
          "invalid input");
      ctd[s[i - 1]] = CTD_3_DANGLE;
    }
    if (c == '5') {
      verify(i + 1 < N && ctd[i + 1] == CTD_NA, "invalid input");
      ctd[i + 1] = CTD_5_DANGLE;
    }
    // Opening a branch.
    if (c == 'p' || c == 'n' || c == '[') {
      if (!stk.empty())
        ++stk.top();  // One more branch here.
      else if (ctd[i] == CTD_NA && c == '[')
        ctd[i] = CTD_UNUSED;
      stk.push(0);  // Number of branches for this part.
    }
    // Closing a branch.
    if (c == 'P' || c == 'N' || c == ']') {
      // Set explicitly unused CTDs for any child branches if this is a multiloop.
      // Also this branch as an outer loop.
      if (stk.top() >= 2) {
        if (ctd[i] == CTD_NA) ctd[i] = CTD_UNUSED;
        for (int j = s[i] + 1; j < i; ++j) {
          if (s[j] == -1) continue;
          if (ctd[j] == CTD_NA) ctd[j] = CTD_UNUSED;
          j = s[j];
        }
      }
      stk.pop();
    }
    if (allowed_characters.find(c) == std::string::npos) error("invalid input '{}'", ctd_str[i]);
  }
  // Add in the terminal mismatches.
  for (int i = 0; i < N; ++i) {
    const char c = ctd_str[i];
    if (c == 'm') {
      const bool right_branch_viable =
          i + 1 < N && s[i + 1] != -1 && s[i + 1] + 1 < N && ctd_str[s[i + 1] + 1] == 'M';
      verify(right_branch_viable && ctd[i + 1] != CTD_NA, "invalid input");
      if (ctd[i + 1] == CTD_UNUSED) ctd[i + 1] = CTD_MISMATCH;
    }
  }
  return {std::move(r), std::move(s), std::move(ctd)};
}

}  // namespace mrna
