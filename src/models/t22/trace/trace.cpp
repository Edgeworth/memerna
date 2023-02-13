// Copyright 2023 Eliot Courtney.
#include "models/t22/trace/trace.h"

#include <algorithm>
#include <compare>
#include <memory>
#include <stack>
#include <variant>

#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/t04/mfe/dp.h"
#include "models/t22/mfe/mfe.h"

namespace mrna::md::t22 {

namespace {

using t04::DP_P;
using t04::DP_U;
using t04::DP_U2;
using t04::DP_U_GU;
using t04::DP_U_RC;
using t04::DP_U_WC;
using t04::EXT;
using t04::EXT_GU;
using t04::EXT_RC;
using t04::EXT_WC;

struct TracebackInternal {
  const Primary& r;
  const Model& em;
  int N;
  const t04::DpArray& dp;
  const t04::ExtArray& ext;
  const Array2D<Energy>& nostack;
  const Array3D<Energy>& penult;
  TraceResult res;
  std::stack<Index> q;

  TracebackInternal(const Primary& r_, const Model::Ptr& em_, const DpState& state_)
      : r(r_), em(*em_), N(static_cast<int>(r_.size())), dp(state_.t04.dp), ext(state_.t04.ext),
        nostack(state_.nostack), penult(state_.penult), res((Secondary(N)), Ctds(N)) {}

  bool ComputeExt(int st, int en, int a) {
    // Case: No pair starting here
    if (a == EXT && st + 1 < N && ext[st + 1][EXT] == ext[st][EXT]) {
      q.emplace(t04::Index(st + 1, -1, EXT));
      return true;
    }
    for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const auto base00 = dp[st][en][DP_P] + em.AuGuPenalty(stb, enb);
      const auto base01 = dp[st][en - 1][DP_P] + em.AuGuPenalty(stb, en1b);
      const auto base10 = dp[st + 1][en][DP_P] + em.AuGuPenalty(st1b, enb);
      const auto base11 = dp[st + 1][en - 1][DP_P] + em.AuGuPenalty(st1b, en1b);

      // (   )<.( * ). > Right coax backward
      if (a == EXT_RC) {
        // Don't set CTDs here since they will have already been set.
        if (base11 + em.MismatchCoaxial(en1b, enb, stb, st1b) + ext[en + 1][EXT] ==
            ext[st][EXT_RC]) {
          q.emplace(t04::Index(st + 1, en - 1, DP_P));
          q.emplace(t04::Index(en + 1, -1, EXT));
          return true;
        }
        continue;
      }

      // (   )<   >
      auto val = base00 + ext[en + 1][EXT];
      if (val == ext[st][a] && (a != EXT_WC || IsWcPair(stb, enb)) &&
          (a != EXT_GU || IsGuPair(stb, enb))) {
        // EXT_WC and EXT_GU will have already had their ctds set.
        if (a == EXT) res.ctd[st] = CTD_UNUSED;
        q.emplace(t04::Index(st, en, DP_P));
        q.emplace(t04::Index(en + 1, -1, EXT));
        return true;
      }

      // Only look at EXT from here on.
      if (a != EXT) continue;

      // (   )3<   > 3'
      if (base01 + em.dangle3[en1b][enb][stb] + ext[en + 1][EXT] == ext[st][EXT]) {
        res.ctd[st] = CTD_3_DANGLE;
        q.emplace(t04::Index(st, en - 1, DP_P));
        q.emplace(t04::Index(en + 1, -1, EXT));
        return true;
      }
      // 5(   )<   > 5'
      if (base10 + em.dangle5[enb][stb][st1b] + ext[en + 1][EXT] == ext[st][EXT]) {
        res.ctd[st + 1] = CTD_5_DANGLE;
        q.emplace(t04::Index(st + 1, en, DP_P));
        q.emplace(t04::Index(en + 1, -1, EXT));
        return true;
      }
      // .(   ).<   > Terminal mismatch
      if (base11 + em.terminal[en1b][enb][stb][st1b] + ext[en + 1][EXT] == ext[st][EXT]) {
        res.ctd[st + 1] = CTD_MISMATCH;
        q.emplace(t04::Index(st + 1, en - 1, DP_P));
        q.emplace(t04::Index(en + 1, -1, EXT));
        return true;
      }

      if (en < N - 1) {
        // .(   ).<(   ) > Left coax  x
        val = base11 + em.MismatchCoaxial(en1b, enb, stb, st1b);
        if (val + ext[en + 1][EXT_WC] == ext[st][EXT]) {
          res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
          res.ctd[en + 1] = CTD_LCOAX_WITH_PREV;
          q.emplace(t04::Index(st + 1, en - 1, DP_P));
          q.emplace(t04::Index(en + 1, -1, EXT_WC));
          return true;
        }
        if (val + ext[en + 1][EXT_GU] == ext[st][EXT]) {
          res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
          res.ctd[en + 1] = CTD_LCOAX_WITH_PREV;
          q.emplace(t04::Index(st + 1, en - 1, DP_P));
          q.emplace(t04::Index(en + 1, -1, EXT_GU));
          return true;
        }

        // (   )<.(   ). > Right coax forward
        if (base00 + ext[en + 1][EXT_RC] == ext[st][EXT]) {
          res.ctd[st] = CTD_RC_WITH_NEXT;
          res.ctd[en + 2] = CTD_RC_WITH_PREV;
          q.emplace(t04::Index(st, en, DP_P));
          q.emplace(t04::Index(en + 1, -1, EXT_RC));
          return true;
        }

        // (   )(<   ) > Flush coax
        if (base01 + em.stack[en1b][enb][WcPair(enb)][stb] + ext[en][EXT_WC] == ext[st][EXT]) {
          res.ctd[st] = CTD_FCOAX_WITH_NEXT;
          res.ctd[en] = CTD_FCOAX_WITH_PREV;
          q.emplace(t04::Index(st, en - 1, DP_P));
          q.emplace(t04::Index(en, -1, EXT_WC));
          return true;
        }
        if (IsGu(enb) &&
            base01 + em.stack[en1b][enb][GuPair(enb)][stb] + ext[en][EXT_GU] == ext[st][EXT]) {
          res.ctd[st] = CTD_FCOAX_WITH_NEXT;
          res.ctd[en] = CTD_FCOAX_WITH_PREV;
          q.emplace(t04::Index(st, en - 1, DP_P));
          q.emplace(t04::Index(en, -1, EXT_GU));
          return true;
        }
      }
    }
    return false;
  }

  bool ComputePairedOrNoStack(int st, int en, bool is_nostack) {
    const auto stb = r[st];
    const auto st1b = r[st + 1];
    const auto st2b = r[st + 2];
    const auto enb = r[en];
    const auto en1b = r[en - 1];
    const auto en2b = r[en - 2];

    // It's paired, so add it to the folding.
    res.s[st] = en;
    res.s[en] = st;

    if (!is_nostack) {
      const int max_stack = en - st - HAIRPIN_MIN_SZ + 1;
      const Energy bulge_left = em.Bulge(r, st, en, st + 2, en - 1);
      const Energy bulge_right = em.Bulge(r, st, en, st + 1, en - 2);

      for (int length = 2; 2 * length <= max_stack; ++length) {
        auto none = em.stack[r[st]][r[st + 1]][r[en - 1]][r[en]] +
            em.penultimate_stack[en1b][enb][stb][st1b];
        if (length == 2 &&
            none + nostack[st + 1][en - 1] +
                    em.penultimate_stack[r[st]][r[st + 1]][r[en - 1]][r[en]] ==
                dp[st][en][DP_P]) {
          q.emplace(NoStackIndex(st + 1, en - 1));
          return true;
        }
        if (none + penult[st + 1][en - 1][length - 1] == dp[st][en][DP_P]) {
          q.emplace(PenultimateIndex(st + 1, en - 1, length - 1));
          return true;
        }

        auto left = bulge_left + em.penultimate_stack[en1b][enb][stb][st2b];
        if (length == 2 &&
            left + nostack[st + 2][en - 1] +
                    em.penultimate_stack[r[st]][r[st + 2]][r[en - 1]][r[en]] ==
                dp[st][en][DP_P]) {
          q.emplace(NoStackIndex(st + 2, en - 1));
          return true;
        }
        if (left + penult[st + 2][en - 1][length - 1] == dp[st][en][DP_P]) {
          q.emplace(PenultimateIndex(st + 2, en - 1, length - 1));
          return true;
        }

        auto right = bulge_right + em.penultimate_stack[en2b][enb][stb][st1b];
        if (length == 2 &&
            right + nostack[st + 1][en - 2] +
                    em.penultimate_stack[r[st]][r[st + 1]][r[en - 2]][r[en]] ==
                dp[st][en][DP_P]) {
          q.emplace(NoStackIndex(st + 1, en - 2));
          return true;
        }
        if (right + penult[st + 1][en - 2][length - 1] == dp[st][en][DP_P]) {
          q.emplace(PenultimateIndex(st + 1, en - 2, length - 1));
          return true;
        }
      }
    }

    const auto target = is_nostack ? nostack[st][en] : dp[st][en][DP_P];

    // Following largely matches the above DP so look up there for comments.
    const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
    for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
      for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
        // Try all internal loops. We don't check stacks or 1 nuc bulge loops.
        if (dp[ist][ien][DP_P] < CAP_E && ist - st + en - ien > 3) {
          const auto val = em.TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P];
          if (val == target) {
            q.emplace(t04::Index(ist, ien, DP_P));
            return true;
          }
        }
      }
    }

    const auto base_branch_cost =
        em.AuGuPenalty(stb, enb) + em.multiloop_hack_a + em.multiloop_hack_b;
    // (<   ><    >)
    if (base_branch_cost + dp[st + 1][en - 1][DP_U2] == target) {
      res.ctd[en] = CTD_UNUSED;
      q.emplace(t04::Index(st + 1, en - 1, DP_U2));
      return true;
    }
    // (3<   ><   >) 3'
    if (base_branch_cost + dp[st + 2][en - 1][DP_U2] + em.dangle3[stb][st1b][enb] == target) {
      res.ctd[en] = CTD_3_DANGLE;
      q.emplace(t04::Index(st + 2, en - 1, DP_U2));
      return true;
    }
    // (<   ><   >5) 5'
    if (base_branch_cost + dp[st + 1][en - 2][DP_U2] + em.dangle5[stb][en1b][enb] == target) {
      res.ctd[en] = CTD_5_DANGLE;
      q.emplace(t04::Index(st + 1, en - 2, DP_U2));
      return true;
    }
    // (.<   ><   >.) Terminal mismatch
    if (base_branch_cost + dp[st + 2][en - 2][DP_U2] + em.terminal[stb][st1b][en1b][enb] ==
        target) {
      res.ctd[en] = CTD_MISMATCH;
      q.emplace(t04::Index(st + 2, en - 2, DP_U2));
      return true;
    }

    for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
      const Base pl1b = r[piv - 1];
      const Base plb = r[piv];
      const Base prb = r[piv + 1];
      const Base pr1b = r[piv + 2];

      // (.(   )   .) Left outer coax - P
      const auto outer_coax = em.MismatchCoaxial(stb, st1b, en1b, enb);
      if (base_branch_cost + dp[st + 2][piv][DP_P] + em.multiloop_hack_b +
              em.AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax ==
          target) {
        res.ctd[en] = CTD_LCOAX_WITH_NEXT;
        res.ctd[st + 2] = CTD_LCOAX_WITH_PREV;
        q.emplace(t04::Index(st + 2, piv, DP_P));
        q.emplace(t04::Index(piv + 1, en - 2, DP_U));
        return true;
      }
      // (.   (   ).) Right outer coax
      if (base_branch_cost + dp[st + 2][piv][DP_U] + em.multiloop_hack_b +
              em.AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax ==
          target) {
        res.ctd[en] = CTD_RC_WITH_PREV;
        res.ctd[piv + 1] = CTD_RC_WITH_NEXT;
        q.emplace(t04::Index(st + 2, piv, DP_U));
        q.emplace(t04::Index(piv + 1, en - 2, DP_P));
        return true;
      }

      // (.(   ).   ) Left inner coax
      if (base_branch_cost + dp[st + 2][piv - 1][DP_P] + em.multiloop_hack_b +
              em.AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
              em.MismatchCoaxial(pl1b, plb, st1b, st2b) ==
          target) {
        res.ctd[en] = CTD_RC_WITH_NEXT;
        res.ctd[st + 2] = CTD_RC_WITH_PREV;
        q.emplace(t04::Index(st + 2, piv - 1, DP_P));
        q.emplace(t04::Index(piv + 1, en - 1, DP_U));
        return true;
      }
      // (   .(   ).) Right inner coax
      if (base_branch_cost + dp[st + 1][piv][DP_U] + em.multiloop_hack_b +
              em.AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
              em.MismatchCoaxial(en2b, en1b, prb, pr1b) ==
          target) {
        res.ctd[en] = CTD_LCOAX_WITH_PREV;
        res.ctd[piv + 2] = CTD_LCOAX_WITH_NEXT;
        q.emplace(t04::Index(st + 1, piv, DP_U));
        q.emplace(t04::Index(piv + 2, en - 2, DP_P));
        return true;
      }

      // ((   )   ) Left flush coax
      if (base_branch_cost + dp[st + 1][piv][DP_P] + em.multiloop_hack_b +
              em.AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] +
              em.stack[stb][st1b][plb][enb] ==
          target) {
        res.ctd[en] = CTD_FCOAX_WITH_NEXT;
        res.ctd[st + 1] = CTD_FCOAX_WITH_PREV;
        q.emplace(t04::Index(st + 1, piv, DP_P));
        q.emplace(t04::Index(piv + 1, en - 1, DP_U));
        return true;
      }
      // (   (   )) Right flush coax
      if (base_branch_cost + dp[st + 1][piv][DP_U] + em.multiloop_hack_b +
              em.AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] +
              em.stack[stb][prb][en1b][enb] ==
          target) {
        res.ctd[en] = CTD_FCOAX_WITH_PREV;
        res.ctd[piv + 1] = CTD_FCOAX_WITH_NEXT;
        q.emplace(t04::Index(st + 1, piv, DP_U));
        q.emplace(t04::Index(piv + 1, en - 1, DP_P));
        return true;
      }
    }

    return false;
  }

  bool ComputeUnpaired(int st, int en, int a) {
    const auto stb = r[st];
    const auto st1b = r[st + 1];

    // Deal with the rest of the cases:
    // Left unpaired. Either DP_U or DP_U2.
    if (st + 1 < en && (a == DP_U || a == DP_U2) && dp[st + 1][en][a] == dp[st][en][a]) {
      q.emplace(t04::Index(st + 1, en, a));
      return true;
    }

    // Pair here.
    for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
      //   (   .   )<   (
      // stb pl1b pb   pr1b
      const auto pb = r[piv];
      const auto pl1b = r[piv - 1];
      // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the
      // right.
      const auto base00 = dp[st][piv][DP_P] + em.AuGuPenalty(stb, pb) + em.multiloop_hack_b;
      const auto base01 = dp[st][piv - 1][DP_P] + em.AuGuPenalty(stb, pl1b) + em.multiloop_hack_b;
      const auto base10 = dp[st + 1][piv][DP_P] + em.AuGuPenalty(st1b, pb) + em.multiloop_hack_b;
      const auto base11 =
          dp[st + 1][piv - 1][DP_P] + em.AuGuPenalty(st1b, pl1b) + em.multiloop_hack_b;

      // Min is for either placing another unpaired or leaving it as nothing.
      // If we're at U2, don't allow leaving as nothing.
      auto right_unpaired = dp[piv + 1][en][DP_U];
      if (a != DP_U2) right_unpaired = std::min(right_unpaired, ZERO_E);

      // Check a == U_RC:
      // (   )<.( ** ). > Right coax backward
      if (a == DP_U_RC) {
        if (base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b) + right_unpaired ==
            dp[st][en][DP_U_RC]) {
          // Ctds were already set from the recurrence that called this.
          q.emplace(t04::Index(st + 1, piv - 1, DP_P));
          if (right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
          return true;
        }
        continue;
      }

      // (   )<   > - U, U2, U_WC?, U_GU?
      if (base00 + right_unpaired == dp[st][en][a] && (a != DP_U_WC || IsWcPair(stb, pb)) &&
          (a != DP_U_GU || IsGuPair(stb, pb))) {
        // If U_WC, or U_GU, we were involved in some sort of coaxial stack previously, and
        // were already set.
        if (a != DP_U_WC && a != DP_U_GU) res.ctd[st] = CTD_UNUSED;
        q.emplace(t04::Index(st, piv, DP_P));
        if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
        return true;
      }

      // The rest of the cases are for U and U2.
      if (a != DP_U && a != DP_U2) continue;

      // (   )3<   > 3' - U, U2
      if (base01 + em.dangle3[pl1b][pb][stb] + right_unpaired == dp[st][en][a]) {
        res.ctd[st] = CTD_3_DANGLE;
        q.emplace(t04::Index(st, piv - 1, DP_P));
        if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
        return true;
      }
      // 5(   )<   > 5' - U, U2
      if (base10 + em.dangle5[pb][stb][st1b] + right_unpaired == dp[st][en][a]) {
        res.ctd[st + 1] = CTD_5_DANGLE;
        q.emplace(t04::Index(st + 1, piv, DP_P));
        if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
        return true;
      }
      // .(   ).<   > Terminal mismatch - U, U2
      if (base11 + em.terminal[pl1b][pb][stb][st1b] + right_unpaired == dp[st][en][a]) {
        res.ctd[st + 1] = CTD_MISMATCH;
        q.emplace(t04::Index(st + 1, piv - 1, DP_P));
        if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
        return true;
      }
      // .(   ).<(   ) > Left coax - U, U2
      auto val = base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b);
      if (val + dp[piv + 1][en][DP_U_WC] == dp[st][en][a]) {
        res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
        res.ctd[piv + 1] = CTD_LCOAX_WITH_PREV;
        q.emplace(t04::Index(st + 1, piv - 1, DP_P));
        q.emplace(t04::Index(piv + 1, en, DP_U_WC));
        return true;
      }
      if (val + dp[piv + 1][en][DP_U_GU] == dp[st][en][a]) {
        res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
        res.ctd[piv + 1] = CTD_LCOAX_WITH_PREV;
        q.emplace(t04::Index(st + 1, piv - 1, DP_P));
        q.emplace(t04::Index(piv + 1, en, DP_U_GU));
        return true;
      }

      // (   )<.(   ). > Right coax forward - U, U2
      if (base00 + dp[piv + 1][en][DP_U_RC] == dp[st][en][a]) {
        res.ctd[st] = CTD_RC_WITH_NEXT;
        res.ctd[piv + 2] = CTD_RC_WITH_PREV;
        q.emplace(t04::Index(st, piv, DP_P));
        q.emplace(t04::Index(piv + 1, en, DP_U_RC));
        return true;
      }

      // (   )(<   ) > Flush coax - U, U2
      if (base01 + em.stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC] == dp[st][en][a]) {
        res.ctd[st] = CTD_FCOAX_WITH_NEXT;
        res.ctd[piv] = CTD_FCOAX_WITH_PREV;
        q.emplace(t04::Index(st, piv - 1, DP_P));
        q.emplace(t04::Index(piv, en, DP_U_WC));
        return true;
      }
      if ((IsGu(pb)) &&
          base01 + em.stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU] == dp[st][en][a]) {
        res.ctd[st] = CTD_FCOAX_WITH_NEXT;
        res.ctd[piv] = CTD_FCOAX_WITH_PREV;
        q.emplace(t04::Index(st, piv - 1, DP_P));
        q.emplace(t04::Index(piv, en, DP_U_GU));
        return true;
      }
    }
    return false;
  }

  bool ComputePenultimate(int st, int en, int length) {
    // It's paired, so add it to the folding.
    res.s[st] = en;
    res.s[en] = st;

    const auto bulge_left = em.Bulge(r, st, en, st + 2, en - 1);
    const auto bulge_right = em.Bulge(r, st, en, st + 1, en - 2);

    auto none = em.stack[r[st]][r[st + 1]][r[en - 1]][r[en]];
    if (length == 2 &&
        none + nostack[st + 1][en - 1] + em.penultimate_stack[r[st]][r[st + 1]][r[en - 1]][r[en]] ==
            penult[st][en][length]) {
      q.emplace(NoStackIndex(st + 1, en - 1));
      return true;
    }
    if (none + penult[st + 1][en - 1][length - 1] == penult[st][en][length]) {
      q.emplace(PenultimateIndex(st + 1, en - 1, length - 1));
      return true;
    }

    auto left = bulge_left;
    if (length == 2 &&
        left + nostack[st + 2][en - 1] + em.penultimate_stack[r[st]][r[st + 2]][r[en - 1]][r[en]] ==
            penult[st][en][length]) {
      q.emplace(NoStackIndex(st + 2, en - 1));
      return true;
    }
    if (left + penult[st + 2][en - 1][length - 1] == penult[st][en][length]) {
      q.emplace(PenultimateIndex(st + 2, en - 1, length - 1));
      return true;
    }

    auto right = bulge_right;
    if (length == 2 &&
        right + nostack[st + 1][en - 2] +
                em.penultimate_stack[r[st]][r[st + 1]][r[en - 2]][r[en]] ==
            penult[st][en][length]) {
      q.emplace(NoStackIndex(st + 1, en - 2));
      return true;
    }
    if (right + penult[st + 1][en - 2][length - 1] == penult[st][en][length]) {
      q.emplace(PenultimateIndex(st + 1, en - 2, length - 1));
      return true;
    }
    return false;
  }

  TraceResult Compute() {
    verify(em.cfg.lonely_pairs != erg::EnergyCfg::LonelyPairs::OFF,
        "fully disallowing lonely pairs is not supported in this energy model");
    verify(em.cfg.ctd == erg::EnergyCfg::Ctd::ALL,
        "only full CTDs are supported in this energy model");

    q.emplace(t04::Index(0, -1, EXT));
    while (!q.empty()) {
      auto idx_all = q.top();
      q.pop();

      if (std::holds_alternative<t04::Index>(idx_all)) {
        auto idx = std::get<t04::Index>(idx_all);
        int st = idx.st;
        int en = idx.en;
        int a = idx.a;
        if (en == -1) {
          ComputeExt(st, en, a);
        } else {
          // Done with paired. We might not have returned true if this was a hairpin.
          if (a == DP_P) {
            ComputePairedOrNoStack(st, en, /*is_nostack=*/false);
          } else {
            [[maybe_unused]] bool found = ComputeUnpaired(st, en, a);
            assert(found);
          }
        }
      } else if (std::holds_alternative<PenultimateIndex>(idx_all)) {
        auto idx = std::get<PenultimateIndex>(idx_all);
        const int st = idx.st;
        const int en = idx.en;
        const int length = idx.len;
        [[maybe_unused]] bool found = ComputePenultimate(st, en, length);
        assert(found);
      } else {
        auto idx = std::get<NoStackIndex>(idx_all);
        const int st = idx.st;
        const int en = idx.en;
        ComputePairedOrNoStack(st, en, /*is_nostack=*/true);
      }
    }

    return std::move(res);
  }
};

}  // namespace

TraceResult Traceback(const Primary& r, const Model::Ptr& em, const DpState& state) {
  return TracebackInternal(r, em, state).Compute();
}

}  // namespace mrna::md::t22
