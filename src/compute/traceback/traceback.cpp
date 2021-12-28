// Copyright 2016 E.
#include "compute/traceback/traceback.h"

#include <algorithm>
#include <stack>

#include "compute/dp.h"
#include "compute/energy/globals.h"
#include "compute/mfe/mfe.h"

namespace mrna::traceback {

// TODO: get rid of these
using energy::gem;
using mfe::internal::gctd;
using mfe::internal::gdp;
using mfe::internal::genergy;
using mfe::internal::gext;
using mfe::internal::gp;

void ComputeExterior() {
  const int N = static_cast<int>(gr.size());
  // Exterior loop calculation. There can be no paired base on gext[en].
  gext[N][EXT] = 0;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    gext[st][EXT] = gext[st + 1][EXT];
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = gr[st], st1b = gr[st + 1], enb = gr[en], en1b = gr[en - 1];
      const auto base00 = gdp[st][en][DP_P] + gem.AuGuPenalty(stb, enb);
      const auto base01 = gdp[st][en - 1][DP_P] + gem.AuGuPenalty(stb, en1b);
      const auto base10 = gdp[st + 1][en][DP_P] + gem.AuGuPenalty(st1b, enb);
      const auto base11 = gdp[st + 1][en - 1][DP_P] + gem.AuGuPenalty(st1b, en1b);
      energy_t ext = MAX_E;

      // (   )<   >
      auto val = base00 + gext[en + 1][EXT];
      ext = std::min(ext, val);
      if (IsGu(stb, enb))
        gext[st][EXT_GU] = std::min(gext[st][EXT_GU], val);
      else
        gext[st][EXT_WC] = std::min(gext[st][EXT_WC], val);

      // (   )3<   > 3'
      ext = std::min(ext, base01 + gem.dangle3[en1b][enb][stb] + gext[en + 1][EXT]);
      // 5(   )<   > 5'
      ext = std::min(ext, base10 + gem.dangle5[enb][stb][st1b] + gext[en + 1][EXT]);
      // .(   ).<   > Terminal mismatch
      ext = std::min(ext, base11 + gem.terminal[en1b][enb][stb][st1b] + gext[en + 1][EXT]);
      // .(   ).<(   ) > Left coax
      val = base11 + gem.MismatchCoaxial(en1b, enb, stb, st1b);
      ext = std::min(ext, val + gext[en + 1][EXT_GU]);
      ext = std::min(ext, val + gext[en + 1][EXT_WC]);

      // (   )<.(   ). > Right coax forward
      ext = std::min(ext, base00 + gext[en + 1][EXT_RCOAX]);
      // (   )<.( * ). > Right coax backward
      gext[st][EXT_RCOAX] = std::min(gext[st][EXT_RCOAX],
          base11 + gem.MismatchCoaxial(en1b, enb, stb, st1b) + gext[en + 1][EXT]);

      // (   )(<   ) > Flush coax
      ext = std::min(ext, base01 + gem.stack[en1b][enb][enb ^ 3][stb] + gext[en][EXT_WC]);
      if (enb == G || enb == U)
        ext = std::min(ext, base01 + gem.stack[en1b][enb][enb ^ 1][stb] + gext[en][EXT_GU]);

      gext[st][EXT] = std::min(gext[st][EXT], ext);
    }
  }
}

void Traceback() {
  const int N = static_cast<int>(gr.size());
  genergy = gext[0][EXT];
  std::stack<index_t> q;
  q.emplace(0, -1, EXT);
  while (!q.empty()) {
    int st = q.top().st, en = q.top().en, a = q.top().a;
    q.pop();

    if (en == -1) {
      // Case: No pair starting here
      if (a == EXT && st + 1 < N && gext[st + 1][EXT] == gext[st][EXT]) {
        q.emplace(st + 1, -1, EXT);
        goto loopend;
      }
      for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        // .   .   .   (   .   .   .   )   <   >
        //           stb  st1b   en1b  enb   rem
        const auto stb = gr[st], st1b = gr[st + 1], enb = gr[en], en1b = gr[en - 1];
        const auto base00 = gdp[st][en][DP_P] + gem.AuGuPenalty(stb, enb);
        const auto base01 = gdp[st][en - 1][DP_P] + gem.AuGuPenalty(stb, en1b);
        const auto base10 = gdp[st + 1][en][DP_P] + gem.AuGuPenalty(st1b, enb);
        const auto base11 = gdp[st + 1][en - 1][DP_P] + gem.AuGuPenalty(st1b, en1b);

        // (   )<.( * ). > Right coax backward
        if (a == EXT_RCOAX) {
          // Don't set CTDs here since they will have already been set.
          if (base11 + gem.MismatchCoaxial(en1b, enb, stb, st1b) + gext[en + 1][EXT] ==
              gext[st][EXT_RCOAX]) {
            q.emplace(st + 1, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT);
            goto loopend;
          }
          continue;
        }

        // (   )<   >
        auto val = base00 + gext[en + 1][EXT];
        if (val == gext[st][a] && (a != EXT_WC || IsWatsonCrick(stb, enb)) &&
            (a != EXT_GU || IsGu(stb, enb))) {
          // EXT_WC and EXT_GU will have already had their ctds set.
          if (a == EXT) gctd[st] = CTD_UNUSED;
          q.emplace(st, en, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }

        // Only look at EXT from here on.
        if (a != EXT) continue;

        // (   )3<   > 3'
        if (base01 + gem.dangle3[en1b][enb][stb] + gext[en + 1][EXT] == gext[st][EXT]) {
          gctd[st] = CTD_3_DANGLE;
          q.emplace(st, en - 1, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }
        // 5(   )<   > 5'
        if (base10 + gem.dangle5[enb][stb][st1b] + gext[en + 1][EXT] == gext[st][EXT]) {
          gctd[st + 1] = CTD_5_DANGLE;
          q.emplace(st + 1, en, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch
        if (base11 + gem.terminal[en1b][enb][stb][st1b] + gext[en + 1][EXT] == gext[st][EXT]) {
          gctd[st + 1] = CTD_MISMATCH;
          q.emplace(st + 1, en - 1, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }

        if (en < N - 1) {
          // .(   ).<(   ) > Left coax  x
          val = base11 + gem.MismatchCoaxial(en1b, enb, stb, st1b);
          if (val + gext[en + 1][EXT_WC] == gext[st][EXT]) {
            gctd[st + 1] = CTD_LCOAX_WITH_NEXT;
            gctd[en + 1] = CTD_LCOAX_WITH_PREV;
            q.emplace(st + 1, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT_WC);
            goto loopend;
          }
          if (val + gext[en + 1][EXT_GU] == gext[st][EXT]) {
            gctd[st + 1] = CTD_LCOAX_WITH_NEXT;
            gctd[en + 1] = CTD_LCOAX_WITH_PREV;
            q.emplace(st + 1, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT_GU);
            goto loopend;
          }

          // (   )<.(   ). > Right coax forward
          if (base00 + gext[en + 1][EXT_RCOAX] == gext[st][EXT]) {
            gctd[st] = CTD_RCOAX_WITH_NEXT;
            gctd[en + 2] = CTD_RCOAX_WITH_PREV;
            q.emplace(st, en, DP_P);
            q.emplace(en + 1, -1, EXT_RCOAX);
            goto loopend;
          }

          // (   )(<   ) > Flush coax
          if (base01 + gem.stack[en1b][enb][enb ^ 3][stb] + gext[en][EXT_WC] == gext[st][EXT]) {
            gctd[st] = CTD_FCOAX_WITH_NEXT;
            gctd[en] = CTD_FCOAX_WITH_PREV;
            q.emplace(st, en - 1, DP_P);
            q.emplace(en, -1, EXT_WC);
            goto loopend;
          }
          if ((enb == G || enb == U) &&
              base01 + gem.stack[en1b][enb][enb ^ 1][stb] + gext[en][EXT_GU] == gext[st][EXT]) {
            gctd[st] = CTD_FCOAX_WITH_NEXT;
            gctd[en] = CTD_FCOAX_WITH_PREV;
            q.emplace(st, en - 1, DP_P);
            q.emplace(en, -1, EXT_GU);
            goto loopend;
          }
        }
      }
    } else {
      const auto stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en],
                 en1b = gr[en - 1], en2b = gr[en - 2];
      if (a == DP_P) {
        // It's paired, so add it to the folding.
        gp[st] = en;
        gp[en] = st;

        // Following largely matches the above DP so look up there for comments.
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (gdp[ist][ien][DP_P] < CAP_E) {
              const auto val = gem.TwoLoop(gr, st, en, ist, ien) + gdp[ist][ien][DP_P];
              if (val == gdp[st][en][DP_P]) {
                q.emplace(ist, ien, DP_P);
                goto loopend;
              }
            }
          }
        }

        const auto base_branch_cost =
            gem.AuGuPenalty(stb, enb) + gem.multiloop_hack_a + gem.multiloop_hack_b;
        // (<   ><    >)
        if (base_branch_cost + gdp[st + 1][en - 1][DP_U2] == gdp[st][en][DP_P]) {
          gctd[en] = CTD_UNUSED;
          q.emplace(st + 1, en - 1, DP_U2);
          goto loopend;
        }
        // (3<   ><   >) 3'
        if (base_branch_cost + gdp[st + 2][en - 1][DP_U2] + gem.dangle3[stb][st1b][enb] ==
            gdp[st][en][DP_P]) {
          gctd[en] = CTD_3_DANGLE;
          q.emplace(st + 2, en - 1, DP_U2);
          goto loopend;
        }
        // (<   ><   >5) 5'
        if (base_branch_cost + gdp[st + 1][en - 2][DP_U2] + gem.dangle5[stb][en1b][enb] ==
            gdp[st][en][DP_P]) {
          gctd[en] = CTD_5_DANGLE;
          q.emplace(st + 1, en - 2, DP_U2);
          goto loopend;
        }
        // (.<   ><   >.) Terminal mismatch
        if (base_branch_cost + gdp[st + 2][en - 2][DP_U2] + gem.terminal[stb][st1b][en1b][enb] ==
            gdp[st][en][DP_P]) {
          gctd[en] = CTD_MISMATCH;
          q.emplace(st + 2, en - 2, DP_U2);
          goto loopend;
        }

        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          const base_t pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1], pr1b = gr[piv + 2];

          // (.(   )   .) Left outer coax - P
          const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);
          if (base_branch_cost + gdp[st + 2][piv][DP_P] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(st2b, plb) + gdp[piv + 1][en - 2][DP_U] + outer_coax ==
              gdp[st][en][DP_P]) {
            gctd[en] = CTD_LCOAX_WITH_NEXT;
            gctd[st + 2] = CTD_LCOAX_WITH_PREV;
            q.emplace(st + 2, piv, DP_P);
            q.emplace(piv + 1, en - 2, DP_U);
            goto loopend;
          }
          // (.   (   ).) Right outer coax
          if (base_branch_cost + gdp[st + 2][piv][DP_U] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(prb, en2b) + gdp[piv + 1][en - 2][DP_P] + outer_coax ==
              gdp[st][en][DP_P]) {
            gctd[en] = CTD_RCOAX_WITH_PREV;
            gctd[piv + 1] = CTD_RCOAX_WITH_NEXT;
            q.emplace(st + 2, piv, DP_U);
            q.emplace(piv + 1, en - 2, DP_P);
            goto loopend;
          }

          // (.(   ).   ) Left right coax
          if (base_branch_cost + gdp[st + 2][piv - 1][DP_P] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(st2b, pl1b) + gdp[piv + 1][en - 1][DP_U] +
                  gem.MismatchCoaxial(pl1b, plb, st1b, st2b) ==
              gdp[st][en][DP_P]) {
            gctd[en] = CTD_RCOAX_WITH_NEXT;
            gctd[st + 2] = CTD_RCOAX_WITH_PREV;
            q.emplace(st + 2, piv - 1, DP_P);
            q.emplace(piv + 1, en - 1, DP_U);
            goto loopend;
          }
          // (   .(   ).) Right left coax
          if (base_branch_cost + gdp[st + 1][piv][DP_U] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(pr1b, en2b) + gdp[piv + 2][en - 2][DP_P] +
                  gem.MismatchCoaxial(en2b, en1b, prb, pr1b) ==
              gdp[st][en][DP_P]) {
            gctd[en] = CTD_LCOAX_WITH_PREV;
            gctd[piv + 2] = CTD_LCOAX_WITH_NEXT;
            q.emplace(st + 1, piv, DP_U);
            q.emplace(piv + 2, en - 2, DP_P);
            goto loopend;
          }

          // ((   )   ) Left flush coax
          if (base_branch_cost + gdp[st + 1][piv][DP_P] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(st1b, plb) + gdp[piv + 1][en - 1][DP_U] +
                  gem.stack[stb][st1b][plb][enb] ==
              gdp[st][en][DP_P]) {
            gctd[en] = CTD_FCOAX_WITH_NEXT;
            gctd[st + 1] = CTD_FCOAX_WITH_PREV;
            q.emplace(st + 1, piv, DP_P);
            q.emplace(piv + 1, en - 1, DP_U);
            goto loopend;
          }
          // (   (   )) Right flush coax
          if (base_branch_cost + gdp[st + 1][piv][DP_U] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(prb, en1b) + gdp[piv + 1][en - 1][DP_P] +
                  gem.stack[stb][prb][en1b][enb] ==
              gdp[st][en][DP_P]) {
            gctd[en] = CTD_FCOAX_WITH_PREV;
            gctd[piv + 1] = CTD_FCOAX_WITH_NEXT;
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
      if (st + 1 < en && (a == DP_U || a == DP_U2) && gdp[st + 1][en][a] == gdp[st][en][a]) {
        q.emplace(st + 1, en, a);
        goto loopend;
      }

      // Pair here.
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = gr[piv], pl1b = gr[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const auto base00 = gdp[st][piv][DP_P] + gem.AuGuPenalty(stb, pb) + gem.multiloop_hack_b;
        const auto base01 =
            gdp[st][piv - 1][DP_P] + gem.AuGuPenalty(stb, pl1b) + gem.multiloop_hack_b;
        const auto base10 =
            gdp[st + 1][piv][DP_P] + gem.AuGuPenalty(st1b, pb) + gem.multiloop_hack_b;
        const auto base11 =
            gdp[st + 1][piv - 1][DP_P] + gem.AuGuPenalty(st1b, pl1b) + gem.multiloop_hack_b;

        // Min is for either placing another unpaired or leaving it as nothing.
        // If we're at U2, don't allow leaving as nothing.
        auto right_unpaired = gdp[piv + 1][en][DP_U];
        if (a != DP_U2) right_unpaired = std::min(right_unpaired, 0);

        // Check a == U_RCOAX:
        // (   )<.( ** ). > Right coax backward
        if (a == DP_U_RCOAX) {
          if (base11 + gem.MismatchCoaxial(pl1b, pb, stb, st1b) + right_unpaired ==
              gdp[st][en][DP_U_RCOAX]) {
            // Ctds were already set from the recurrence that called this.
            q.emplace(st + 1, piv - 1, DP_P);
            if (right_unpaired) q.emplace(piv + 1, en, DP_U);
            goto loopend;
          }
          continue;
        }

        // (   )<   > - U, U2, U_WC?, U_GU?
        if (base00 + right_unpaired == gdp[st][en][a] && (a != DP_U_WC || IsWatsonCrick(stb, pb)) &&
            (a != DP_U_GU || IsGu(stb, pb))) {
          // If U_WC, or U_GU, we were involved in some sort of coaxial stack previously, and were
          // already set.
          if (a != DP_U_WC && a != DP_U_GU) gctd[st] = CTD_UNUSED;
          q.emplace(st, piv, DP_P);
          if (a == DP_U2 || right_unpaired) q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }

        // The rest of the cases are for U and U2.
        if (a != DP_U && a != DP_U2) continue;

        // (   )3<   > 3' - U, U2
        if (base01 + gem.dangle3[pl1b][pb][stb] + right_unpaired == gdp[st][en][a]) {
          gctd[st] = CTD_3_DANGLE;
          q.emplace(st, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired) q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // 5(   )<   > 5' - U, U2
        if (base10 + gem.dangle5[pb][stb][st1b] + right_unpaired == gdp[st][en][a]) {
          gctd[st + 1] = CTD_5_DANGLE;
          q.emplace(st + 1, piv, DP_P);
          if (a == DP_U2 || right_unpaired) q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch - U, U2
        if (base11 + gem.terminal[pl1b][pb][stb][st1b] + right_unpaired == gdp[st][en][a]) {
          gctd[st + 1] = CTD_MISMATCH;
          q.emplace(st + 1, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired) q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<(   ) > Left coax - U, U2
        auto val = base11 + gem.MismatchCoaxial(pl1b, pb, stb, st1b);
        if (val + gdp[piv + 1][en][DP_U_WC] == gdp[st][en][a]) {
          gctd[st + 1] = CTD_LCOAX_WITH_NEXT;
          gctd[piv + 1] = CTD_LCOAX_WITH_PREV;
          q.emplace(st + 1, piv - 1, DP_P);
          q.emplace(piv + 1, en, DP_U_WC);
          goto loopend;
        }
        if (val + gdp[piv + 1][en][DP_U_GU] == gdp[st][en][a]) {
          gctd[st + 1] = CTD_LCOAX_WITH_NEXT;
          gctd[piv + 1] = CTD_LCOAX_WITH_PREV;
          q.emplace(st + 1, piv - 1, DP_P);
          q.emplace(piv + 1, en, DP_U_GU);
          goto loopend;
        }

        // (   )<.(   ). > Right coax forward - U, U2
        if (base00 + gdp[piv + 1][en][DP_U_RCOAX] == gdp[st][en][a]) {
          gctd[st] = CTD_RCOAX_WITH_NEXT;
          gctd[piv + 2] = CTD_RCOAX_WITH_PREV;
          q.emplace(st, piv, DP_P);
          q.emplace(piv + 1, en, DP_U_RCOAX);
          goto loopend;
        }

        // (   )(<   ) > Flush coax - U, U2
        if (base01 + gem.stack[pl1b][pb][pb ^ 3][stb] + gdp[piv][en][DP_U_WC] == gdp[st][en][a]) {
          gctd[st] = CTD_FCOAX_WITH_NEXT;
          gctd[piv] = CTD_FCOAX_WITH_PREV;
          q.emplace(st, piv - 1, DP_P);
          q.emplace(piv, en, DP_U_WC);
          goto loopend;
        }
        if ((pb == G || pb == U) &&
            base01 + gem.stack[pl1b][pb][pb ^ 1][stb] + gdp[piv][en][DP_U_GU] == gdp[st][en][a]) {
          gctd[st] = CTD_FCOAX_WITH_NEXT;
          gctd[piv] = CTD_FCOAX_WITH_PREV;
          q.emplace(st, piv - 1, DP_P);
          q.emplace(piv, en, DP_U_GU);
          goto loopend;
        }
      }
    }

  loopend : {};
  }
}

}  // namespace mrna::traceback