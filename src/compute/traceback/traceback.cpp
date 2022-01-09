// Copyright 2016 E.
#include "compute/traceback/traceback.h"

#include <algorithm>
#include <stack>

#include "compute/dp.h"
#include "compute/mfe/mfe.h"

namespace mrna::traceback {

std::tuple<std::vector<int>, std::vector<Ctd>> Traceback(
    const Primary& r, const energy::EnergyModel& em, const DpArray& dp, const ExtArray& ext) {
  const int N = static_cast<int>(r.size());
  std::vector<int> p(N, -1);
  std::vector<Ctd> ctd(N, CTD_NA);

  std::stack<Index> q;
  q.emplace(0, -1, EXT);
  while (!q.empty()) {
    int st = q.top().st, en = q.top().en, a = q.top().a;
    q.pop();

    if (en == -1) {
      // Case: No pair starting here
      if (a == EXT && st + 1 < N && ext[st + 1][EXT] == ext[st][EXT]) {
        q.emplace(st + 1, -1, EXT);
        goto loopend;
      }
      for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        // .   .   .   (   .   .   .   )   <   >
        //           stb  st1b   en1b  enb   rem
        const auto stb = r[st], st1b = r[st + 1], enb = r[en], en1b = r[en - 1];
        const auto base00 = dp[st][en][DP_P] + em.AuGuPenalty(stb, enb);
        const auto base01 = dp[st][en - 1][DP_P] + em.AuGuPenalty(stb, en1b);
        const auto base10 = dp[st + 1][en][DP_P] + em.AuGuPenalty(st1b, enb);
        const auto base11 = dp[st + 1][en - 1][DP_P] + em.AuGuPenalty(st1b, en1b);

        // (   )<.( * ). > Right coax backward
        if (a == EXT_RCOAX) {
          // Don't set CTDs here since they will have already been set.
          if (base11 + em.MismatchCoaxial(en1b, enb, stb, st1b) + ext[en + 1][EXT] ==
              ext[st][EXT_RCOAX]) {
            q.emplace(st + 1, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT);
            goto loopend;
          }
          continue;
        }

        // (   )<   >
        auto val = base00 + ext[en + 1][EXT];
        if (val == ext[st][a] && (a != EXT_WC || IsWatsonCrick(stb, enb)) &&
            (a != EXT_GU || IsGu(stb, enb))) {
          // EXT_WC and EXT_GU will have already had their ctds set.
          if (a == EXT) ctd[st] = CTD_UNUSED;
          q.emplace(st, en, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }

        // Only look at EXT from here on.
        if (a != EXT) continue;

        // (   )3<   > 3'
        if (base01 + em.dangle3[en1b][enb][stb] + ext[en + 1][EXT] == ext[st][EXT]) {
          ctd[st] = CTD_3_DANGLE;
          q.emplace(st, en - 1, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }
        // 5(   )<   > 5'
        if (base10 + em.dangle5[enb][stb][st1b] + ext[en + 1][EXT] == ext[st][EXT]) {
          ctd[st + 1] = CTD_5_DANGLE;
          q.emplace(st + 1, en, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch
        if (base11 + em.terminal[en1b][enb][stb][st1b] + ext[en + 1][EXT] == ext[st][EXT]) {
          ctd[st + 1] = CTD_MISMATCH;
          q.emplace(st + 1, en - 1, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }

        if (en < N - 1) {
          // .(   ).<(   ) > Left coax  x
          val = base11 + em.MismatchCoaxial(en1b, enb, stb, st1b);
          if (val + ext[en + 1][EXT_WC] == ext[st][EXT]) {
            ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
            ctd[en + 1] = CTD_LCOAX_WITH_PREV;
            q.emplace(st + 1, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT_WC);
            goto loopend;
          }
          if (val + ext[en + 1][EXT_GU] == ext[st][EXT]) {
            ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
            ctd[en + 1] = CTD_LCOAX_WITH_PREV;
            q.emplace(st + 1, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT_GU);
            goto loopend;
          }

          // (   )<.(   ). > Right coax forward
          if (base00 + ext[en + 1][EXT_RCOAX] == ext[st][EXT]) {
            ctd[st] = CTD_RCOAX_WITH_NEXT;
            ctd[en + 2] = CTD_RCOAX_WITH_PREV;
            q.emplace(st, en, DP_P);
            q.emplace(en + 1, -1, EXT_RCOAX);
            goto loopend;
          }

          // (   )(<   ) > Flush coax
          if (base01 + em.stack[en1b][enb][enb ^ 3][stb] + ext[en][EXT_WC] == ext[st][EXT]) {
            ctd[st] = CTD_FCOAX_WITH_NEXT;
            ctd[en] = CTD_FCOAX_WITH_PREV;
            q.emplace(st, en - 1, DP_P);
            q.emplace(en, -1, EXT_WC);
            goto loopend;
          }
          if ((enb == G || enb == U) &&
              base01 + em.stack[en1b][enb][enb ^ 1][stb] + ext[en][EXT_GU] == ext[st][EXT]) {
            ctd[st] = CTD_FCOAX_WITH_NEXT;
            ctd[en] = CTD_FCOAX_WITH_PREV;
            q.emplace(st, en - 1, DP_P);
            q.emplace(en, -1, EXT_GU);
            goto loopend;
          }
        }
      }
    } else {
      const auto stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1],
                 en2b = r[en - 2];
      if (a == DP_P) {
        // It's paired, so add it to the folding.
        p[st] = en;
        p[en] = st;

        // Following largely matches the above DP so look up there for comments.
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (dp[ist][ien][DP_P] < CAP_E) {
              const auto val = em.TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P];
              if (val == dp[st][en][DP_P]) {
                q.emplace(ist, ien, DP_P);
                goto loopend;
              }
            }
          }
        }

        const auto base_branch_cost =
            em.AuGuPenalty(stb, enb) + em.multiloop_hack_a + em.multiloop_hack_b;
        // (<   ><    >)
        if (base_branch_cost + dp[st + 1][en - 1][DP_U2] == dp[st][en][DP_P]) {
          ctd[en] = CTD_UNUSED;
          q.emplace(st + 1, en - 1, DP_U2);
          goto loopend;
        }
        // (3<   ><   >) 3'
        if (base_branch_cost + dp[st + 2][en - 1][DP_U2] + em.dangle3[stb][st1b][enb] ==
            dp[st][en][DP_P]) {
          ctd[en] = CTD_3_DANGLE;
          q.emplace(st + 2, en - 1, DP_U2);
          goto loopend;
        }
        // (<   ><   >5) 5'
        if (base_branch_cost + dp[st + 1][en - 2][DP_U2] + em.dangle5[stb][en1b][enb] ==
            dp[st][en][DP_P]) {
          ctd[en] = CTD_5_DANGLE;
          q.emplace(st + 1, en - 2, DP_U2);
          goto loopend;
        }
        // (.<   ><   >.) Terminal mismatch
        if (base_branch_cost + dp[st + 2][en - 2][DP_U2] + em.terminal[stb][st1b][en1b][enb] ==
            dp[st][en][DP_P]) {
          ctd[en] = CTD_MISMATCH;
          q.emplace(st + 2, en - 2, DP_U2);
          goto loopend;
        }

        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          const Base pl1b = r[piv - 1], plb = r[piv], prb = r[piv + 1], pr1b = r[piv + 2];

          // (.(   )   .) Left outer coax - P
          const auto outer_coax = em.MismatchCoaxial(stb, st1b, en1b, enb);
          if (base_branch_cost + dp[st + 2][piv][DP_P] + em.multiloop_hack_b +
                  em.AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax ==
              dp[st][en][DP_P]) {
            ctd[en] = CTD_LCOAX_WITH_NEXT;
            ctd[st + 2] = CTD_LCOAX_WITH_PREV;
            q.emplace(st + 2, piv, DP_P);
            q.emplace(piv + 1, en - 2, DP_U);
            goto loopend;
          }
          // (.   (   ).) Right outer coax
          if (base_branch_cost + dp[st + 2][piv][DP_U] + em.multiloop_hack_b +
                  em.AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax ==
              dp[st][en][DP_P]) {
            ctd[en] = CTD_RCOAX_WITH_PREV;
            ctd[piv + 1] = CTD_RCOAX_WITH_NEXT;
            q.emplace(st + 2, piv, DP_U);
            q.emplace(piv + 1, en - 2, DP_P);
            goto loopend;
          }

          // (.(   ).   ) Left right coax
          if (base_branch_cost + dp[st + 2][piv - 1][DP_P] + em.multiloop_hack_b +
                  em.AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
                  em.MismatchCoaxial(pl1b, plb, st1b, st2b) ==
              dp[st][en][DP_P]) {
            ctd[en] = CTD_RCOAX_WITH_NEXT;
            ctd[st + 2] = CTD_RCOAX_WITH_PREV;
            q.emplace(st + 2, piv - 1, DP_P);
            q.emplace(piv + 1, en - 1, DP_U);
            goto loopend;
          }
          // (   .(   ).) Right left coax
          if (base_branch_cost + dp[st + 1][piv][DP_U] + em.multiloop_hack_b +
                  em.AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
                  em.MismatchCoaxial(en2b, en1b, prb, pr1b) ==
              dp[st][en][DP_P]) {
            ctd[en] = CTD_LCOAX_WITH_PREV;
            ctd[piv + 2] = CTD_LCOAX_WITH_NEXT;
            q.emplace(st + 1, piv, DP_U);
            q.emplace(piv + 2, en - 2, DP_P);
            goto loopend;
          }

          // ((   )   ) Left flush coax
          if (base_branch_cost + dp[st + 1][piv][DP_P] + em.multiloop_hack_b +
                  em.AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] +
                  em.stack[stb][st1b][plb][enb] ==
              dp[st][en][DP_P]) {
            ctd[en] = CTD_FCOAX_WITH_NEXT;
            ctd[st + 1] = CTD_FCOAX_WITH_PREV;
            q.emplace(st + 1, piv, DP_P);
            q.emplace(piv + 1, en - 1, DP_U);
            goto loopend;
          }
          // (   (   )) Right flush coax
          if (base_branch_cost + dp[st + 1][piv][DP_U] + em.multiloop_hack_b +
                  em.AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] +
                  em.stack[stb][prb][en1b][enb] ==
              dp[st][en][DP_P]) {
            ctd[en] = CTD_FCOAX_WITH_PREV;
            ctd[piv + 1] = CTD_FCOAX_WITH_NEXT;
            q.emplace(st + 1, piv, DP_U);
            q.emplace(piv + 1, en - 1, DP_P);
            goto loopend;
          }
        }

        // Done with paired. We might not have jumped to loopend if this was a hairpin.
        continue;
      }

      // Deal with the rest of the cases:
      // Left unpaired. Either DP_U or DP_U2.
      if (st + 1 < en && (a == DP_U || a == DP_U2) && dp[st + 1][en][a] == dp[st][en][a]) {
        q.emplace(st + 1, en, a);
        goto loopend;
      }

      // Pair here.
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = r[piv], pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const auto base00 = dp[st][piv][DP_P] + em.AuGuPenalty(stb, pb) + em.multiloop_hack_b;
        const auto base01 = dp[st][piv - 1][DP_P] + em.AuGuPenalty(stb, pl1b) + em.multiloop_hack_b;
        const auto base10 = dp[st + 1][piv][DP_P] + em.AuGuPenalty(st1b, pb) + em.multiloop_hack_b;
        const auto base11 =
            dp[st + 1][piv - 1][DP_P] + em.AuGuPenalty(st1b, pl1b) + em.multiloop_hack_b;

        // Min is for either placing another unpaired or leaving it as nothing.
        // If we're at U2, don't allow leaving as nothing.
        auto right_unpaired = dp[piv + 1][en][DP_U];
        if (a != DP_U2) right_unpaired = std::min(right_unpaired, 0);

        // Check a == U_RCOAX:
        // (   )<.( ** ). > Right coax backward
        if (a == DP_U_RCOAX) {
          if (base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b) + right_unpaired ==
              dp[st][en][DP_U_RCOAX]) {
            // Ctds were already set from the recurrence that called this.
            q.emplace(st + 1, piv - 1, DP_P);
            if (right_unpaired) q.emplace(piv + 1, en, DP_U);
            goto loopend;
          }
          continue;
        }

        // (   )<   > - U, U2, U_WC?, U_GU?
        if (base00 + right_unpaired == dp[st][en][a] && (a != DP_U_WC || IsWatsonCrick(stb, pb)) &&
            (a != DP_U_GU || IsGu(stb, pb))) {
          // If U_WC, or U_GU, we were involved in some sort of coaxial stack previously, and were
          // already set.
          if (a != DP_U_WC && a != DP_U_GU) ctd[st] = CTD_UNUSED;
          q.emplace(st, piv, DP_P);
          if (a == DP_U2 || right_unpaired) q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }

        // The rest of the cases are for U and U2.
        if (a != DP_U && a != DP_U2) continue;

        // (   )3<   > 3' - U, U2
        if (base01 + em.dangle3[pl1b][pb][stb] + right_unpaired == dp[st][en][a]) {
          ctd[st] = CTD_3_DANGLE;
          q.emplace(st, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired) q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // 5(   )<   > 5' - U, U2
        if (base10 + em.dangle5[pb][stb][st1b] + right_unpaired == dp[st][en][a]) {
          ctd[st + 1] = CTD_5_DANGLE;
          q.emplace(st + 1, piv, DP_P);
          if (a == DP_U2 || right_unpaired) q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch - U, U2
        if (base11 + em.terminal[pl1b][pb][stb][st1b] + right_unpaired == dp[st][en][a]) {
          ctd[st + 1] = CTD_MISMATCH;
          q.emplace(st + 1, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired) q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<(   ) > Left coax - U, U2
        auto val = base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b);
        if (val + dp[piv + 1][en][DP_U_WC] == dp[st][en][a]) {
          ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
          ctd[piv + 1] = CTD_LCOAX_WITH_PREV;
          q.emplace(st + 1, piv - 1, DP_P);
          q.emplace(piv + 1, en, DP_U_WC);
          goto loopend;
        }
        if (val + dp[piv + 1][en][DP_U_GU] == dp[st][en][a]) {
          ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
          ctd[piv + 1] = CTD_LCOAX_WITH_PREV;
          q.emplace(st + 1, piv - 1, DP_P);
          q.emplace(piv + 1, en, DP_U_GU);
          goto loopend;
        }

        // (   )<.(   ). > Right coax forward - U, U2
        if (base00 + dp[piv + 1][en][DP_U_RCOAX] == dp[st][en][a]) {
          ctd[st] = CTD_RCOAX_WITH_NEXT;
          ctd[piv + 2] = CTD_RCOAX_WITH_PREV;
          q.emplace(st, piv, DP_P);
          q.emplace(piv + 1, en, DP_U_RCOAX);
          goto loopend;
        }

        // (   )(<   ) > Flush coax - U, U2
        if (base01 + em.stack[pl1b][pb][pb ^ 3][stb] + dp[piv][en][DP_U_WC] == dp[st][en][a]) {
          ctd[st] = CTD_FCOAX_WITH_NEXT;
          ctd[piv] = CTD_FCOAX_WITH_PREV;
          q.emplace(st, piv - 1, DP_P);
          q.emplace(piv, en, DP_U_WC);
          goto loopend;
        }
        if ((pb == G || pb == U) &&
            base01 + em.stack[pl1b][pb][pb ^ 1][stb] + dp[piv][en][DP_U_GU] == dp[st][en][a]) {
          ctd[st] = CTD_FCOAX_WITH_NEXT;
          ctd[piv] = CTD_FCOAX_WITH_PREV;
          q.emplace(st, piv - 1, DP_P);
          q.emplace(piv, en, DP_U_GU);
          goto loopend;
        }
      }
    }

  loopend : {};
  }

  return {std::move(p), std::move(ctd)};
}

}  // namespace mrna::traceback
