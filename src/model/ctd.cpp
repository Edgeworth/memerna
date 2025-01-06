// Copyright 2022 Eliot Courtney.
#include "model/ctd.h"

#include <utility>
#include <vector>

#include "model/branch.h"
#include "model/primary.h"
#include "util/error.h"

namespace mrna {

namespace {

Ctds MakeD2CtdsForSecondary(const Secondary& s) {
  const int N = static_cast<int>(s.size());
  Ctds ctd(N);
  auto counts = GetBranchCounts(s);
  std::vector<bool> in_multiloop;
  in_multiloop.push_back(true);  // Count exterior loop as multiloop.

  for (int i = 0; i < N; ++i) {
    if (s[i] == -1) continue;
    auto cur_ctd = CTD_UNUSED;
    if (i != 0 && s[i] != N - 1) {
      cur_ctd = CTD_MISMATCH;
    } else if (s[i] != N - 1) {
      cur_ctd = CTD_3_DANGLE;
    } else if (i != 0) {
      cur_ctd = CTD_5_DANGLE;
    }

    // Set CTDs when we are in a multiloop for the outer face, or if we are a multiloop for the
    // inner face.
    if (in_multiloop.back()) ctd[i] = cur_ctd;

    if (s[i] > i)
      in_multiloop.push_back(counts[s[i]] > 1);
    else
      in_multiloop.pop_back();
  }
  return ctd;
}

Ctds ParseD2String(const Secondary& s, const std::string& ctd_str) {
  const int N = static_cast<int>(ctd_str.size());
  Ctds ctd = MakeD2CtdsForSecondary(s);
  for (int i = 0; i < N; ++i) {
    const char c = ctd_str[i];
    verify(c == '[' || c == ']' || c == '.', "invalid input");
    if (c == '.') verify(ctd[i] == CTD_NA, "invalid input");
  }
  return ctd;
}

// Breaking ambiguous case for the old format:
// >(   )>   )<<   )(   )>   )<   )
// (m(   )m(   )m(   )m)
Ctds ParseCtdString(const Secondary& s, const std::string& ctd_str) {
  const int N = static_cast<int>(ctd_str.size());
  const std::string allowed_characters = "[]nNpP.mM35";
  Ctds ctd(N);
  std::vector<int> stk;
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
      // If we optimistically set an exterior loop to CTD_UNUSED, we might want to rewrite it
      // here.
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
        ++stk.back();  // One more branch here.
      else if (ctd[i] == CTD_NA && c == '[')
        ctd[i] = CTD_UNUSED;
      stk.push_back(0);  // Number of branches for this part.
    }
    // Closing a branch.
    if (c == 'P' || c == 'N' || c == ']') {
      // Set explicitly unused CTDs for any child branches if this is a multiloop.
      // Also this branch as an outer loop.
      if (stk.back() >= 2) {
        if (ctd[i] == CTD_NA) ctd[i] = CTD_UNUSED;
        for (int j = s[i] + 1; j < i; ++j) {
          if (s[j] == -1) continue;
          if (ctd[j] == CTD_NA) ctd[j] = CTD_UNUSED;
          j = s[j];
        }
      }
      stk.pop_back();
    }
    if (allowed_characters.find(c) == std::string::npos) fatal("invalid input '{}'", ctd_str[i]);
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
  return ctd;
}

}  // namespace

std::ostream& operator<<(std::ostream& str, const Ctd& o) { return str << CtdToName(o); }

const char* CtdToName(Ctd ctd) {
  switch (ctd) {
  case CTD_NA: return "n/a";
  case CTD_UNUSED: return "unused";
  case CTD_3_DANGLE: return "3' dangle";
  case CTD_5_DANGLE: return "5' dangle";
  case CTD_MISMATCH: return "terminal mismatch";
  case CTD_LCOAX_WITH_NEXT: return "left mismatch coax with next";
  case CTD_LCOAX_WITH_PREV: return "left mismatch coax with prev";
  case CTD_RC_WITH_NEXT: return "right mismatch coax with next";
  case CTD_RC_WITH_PREV: return "right mismatch coax with prev";
  case CTD_FCOAX_WITH_NEXT: return "flush coax with next";
  case CTD_FCOAX_WITH_PREV: return "flush coax with prev";
  case CTD_SIZE: break;
  }
  unreachable();
}

std::string Ctds::ToString(const Secondary& s, bool use_d2) const {
  int N = static_cast<int>(data_.size());
  std::string str(s.size(), '.');

  if (use_d2) {
    Ctds ctd = MakeD2CtdsForSecondary(s);
    for (int i = 0; i < N; ++i) {
      if (s[i] == -1) continue;
      const bool closing = s[i] < i;
      if (closing)
        str[i] = ']';
      else
        str[i] = '[';
      verify(data_[i] == ctd[i], "invalid input: {} != {}", data_[i], ctd[i]);
    }
    return str;
  }

  for (int i = 0; i < N; ++i) {
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
    case CTD_SIZE: unreachable();
    }
  }
  return str;
}

bool Ctds::IsCtdString(const std::string& ctd_str) {
  return std::ranges::none_of(ctd_str, [](char c) { return c == '('; });
}

std::tuple<Primary, Secondary, Ctds> ParseSeqCtdString(
    const std::string& prim_str, const std::string& ctd_str, bool use_d2) {
  auto r = Primary::FromSeq(prim_str);
  verify(prim_str.size() == ctd_str.size(), "primary and pairs string need to be the same length");
  const int N = static_cast<int>(r.size());
  std::vector<int> stk;
  Secondary s(N);
  for (int i = 0; i < N; ++i) {
    const char c = ctd_str[i];
    if (c == 'p' || c == 'n' || c == '[') {
      stk.push_back(i);
    } else if (c == 'P' || c == 'N' || c == ']') {
      verify(!stk.empty(), "invalid input");
      s[stk.back()] = i;
      s[i] = stk.back();
      stk.pop_back();
    }
  }
  verify(stk.empty(), "invalid input");
  auto ctd = use_d2 ? ParseD2String(s, ctd_str) : ParseCtdString(s, ctd_str);

  return {std::move(r), std::move(s), std::move(ctd)};
}

}  // namespace mrna
