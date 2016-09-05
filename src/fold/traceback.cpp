#include "fold/fold.h"

namespace memerna {
namespace fold {

using namespace constants;
using namespace energy;
using namespace internal;

#define UPDATE_EXT(a_, na_, value_) \
  do { \
    energy_t macro_upd_value_ = (value_) + exterior[en + 1][(na_)]; \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < exterior[st][(a_)]) \
      exterior[st][(a_)] = macro_upd_value_; \
  } while (0)

void Context::ComputeExterior() {
  // Exterior loop calculation. There can be no paired base on exterior[en].
  exterior[N][EXT] = 0;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    exterior[st][EXT] = exterior[st + 1][EXT];
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = r[st], st1b = r[st + 1], enb = r[en], en1b = r[en - 1];
      const auto base00 = arr[st][en][DP_P] + em.AuGuPenalty(stb, enb);
      const auto base01 = arr[st][en - 1][DP_P] + em.AuGuPenalty(stb, en1b);
      const auto base10 = arr[st + 1][en][DP_P] + em.AuGuPenalty(st1b, enb);
      const auto base11 = arr[st + 1][en - 1][DP_P] + em.AuGuPenalty(st1b, en1b);

      // (   )<   >
      UPDATE_EXT(EXT, EXT, base00);
      if (IsGu(stb, enb))
        UPDATE_EXT(EXT_GU, EXT, base00);
      else
        UPDATE_EXT(EXT_WC, EXT, base00);

      // (   )3<   > 3'
      UPDATE_EXT(EXT, EXT, base01 + em.dangle3[en1b][enb][stb]);
      // 5(   )<   > 5'
      UPDATE_EXT(EXT, EXT, base10 + em.dangle5[enb][stb][st1b]);
      // .(   ).<   > Terminal mismatch
      UPDATE_EXT(EXT, EXT, base11 + em.terminal[en1b][enb][stb][st1b]);
      // .(   ).<(   ) > Left coax  x
      auto val = base11 + em.MismatchCoaxial(en1b, enb, stb, st1b);
      UPDATE_EXT(EXT, EXT_GU, val);
      UPDATE_EXT(EXT, EXT_WC, val);

      // (   ).<(   ). > Right coax forwardconst
      UPDATE_EXT(EXT, EXT_RCOAX, base01);
      // (   ).<( * ). > Right coax backward
      if (st > 0)
        UPDATE_EXT(EXT_RCOAX, EXT, base01 + em.MismatchCoaxial(en1b, enb, r[st - 1], stb));

      if (en < N - 1) {
        // (   )<(   ) > Flush coax
        const auto enrb = r[en + 1];
        UPDATE_EXT(EXT, EXT_WC, base00 + em.stack[enb][enrb][enrb ^ 3][stb]);
        if (enrb == G || enrb == U)
          UPDATE_EXT(EXT, EXT_GU, base00 + em.stack[enb][enrb][enrb ^ 1][stb]);
      }
    }
  }
}

#undef UPDATE_EXT

computed_t Context::Traceback() {
  computed_t computed(r);
  computed.energy = exterior[0][EXT];
  std::stack<index_t> q;
  q.emplace(0, -1, EXT);
  while (!q.empty()) {
    int st = q.top().st, en = q.top().en, a = q.top().a;
    q.pop();

    if (en == -1) {
      // Case: No pair starting here
      if (a == EXT && st + 1 < N && exterior[st + 1][EXT] == exterior[st][EXT]) {
        q.emplace(st + 1, -1, EXT);
        goto loopend;
      }
      for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        // .   .   .   (   .   .   .   )   <   >
        //           stb  st1b   en1b  enb   rem
        const auto stb = r[st], st1b = r[st + 1], enb = r[en], en1b = r[en - 1];
        const auto base00 = arr[st][en][DP_P] + em.AuGuPenalty(stb, enb);
        const auto base01 = arr[st][en - 1][DP_P] + em.AuGuPenalty(stb, en1b);
        const auto base10 = arr[st + 1][en][DP_P] + em.AuGuPenalty(st1b, enb);
        const auto base11 = arr[st + 1][en - 1][DP_P] + em.AuGuPenalty(st1b, en1b);

        // (   ).<( * ). > Right coax backward
        if (a == EXT_RCOAX) {
          // Don't set CTDs here since they will have already been set.
          if (st > 0 && base01 + em.MismatchCoaxial(en1b, enb, r[st - 1], stb) +
              exterior[en + 1][EXT] == exterior[st][EXT_RCOAX]) {
            q.emplace(st, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT);
            goto loopend;
          }
          continue;
        }

        // (   )<   >
        auto val = base00 + exterior[en + 1][EXT];
        if (val == exterior[st][a] && (a != EXT_WC || IsWatsonCrick(stb, enb)) && (a != EXT_GU || IsGu(stb, enb))) {
          // EXT_WC and EXT_GU will have already had their ctds set.
          if (a == EXT)
            computed.base_ctds[st] = CTD_UNUSED;
          q.emplace(st, en, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }

        // Only look at EXT from here on.
        if (a != EXT) continue;

        // (   )3<   > 3'
        if (base01 + em.dangle3[en1b][enb][stb] + exterior[en + 1][EXT] == exterior[st][EXT]) {
          computed.base_ctds[st] = CTD_3_DANGLE;
          q.emplace(st, en - 1, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }
        // 5(   )<   > 5'
        if (base10 + em.dangle5[enb][stb][st1b] + exterior[en + 1][EXT] == exterior[st][EXT]) {
          computed.base_ctds[st + 1] = CTD_5_DANGLE;
          q.emplace(st + 1, en, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch
        if (base11 + em.terminal[en1b][enb][stb][st1b] + exterior[en + 1][EXT] == exterior[st][EXT]) {
          computed.base_ctds[st + 1] = CTD_TERMINAL_MISMATCH;
          q.emplace(st + 1, en - 1, DP_P);
          q.emplace(en + 1, -1, EXT);
          goto loopend;
        }

        if (en < N - 1) {
          // .(   ).<(   ) > Left coax  x
          val = base11 + em.MismatchCoaxial(en1b, enb, stb, st1b);
          if (val + exterior[en + 1][EXT_WC] == exterior[st][EXT]) {
            computed.base_ctds[st + 1] = CTD_LEFT_MISMATCH_COAX_WITH_NEXT;
            computed.base_ctds[en + 1] = CTD_LEFT_MISMATCH_COAX_WITH_PREV;
            q.emplace(st + 1, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT_WC);
            goto loopend;
          }
          if (val + exterior[en + 1][EXT_GU] == exterior[st][EXT]) {
            computed.base_ctds[st + 1] = CTD_LEFT_MISMATCH_COAX_WITH_NEXT;
            computed.base_ctds[en + 1] = CTD_LEFT_MISMATCH_COAX_WITH_PREV;
            q.emplace(st + 1, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT_GU);
            goto loopend;
          }

          // (   ).<(   ). > Right coax forward
          if (base01 + exterior[en + 1][EXT_RCOAX] == exterior[st][EXT]) {
            computed.base_ctds[st] = CTD_RIGHT_MISMATCH_COAX_WITH_NEXT;
            computed.base_ctds[en + 1] = CTD_RIGHT_MISMATCH_COAX_WITH_PREV;
            q.emplace(st, en - 1, DP_P);
            q.emplace(en + 1, -1, EXT_RCOAX);
            goto loopend;
          }

          // (   )<(   ) > Flush coax
          const auto enrb = r[en + 1];
          if (base00 + em.stack[enb][enrb][enrb ^ 3][stb] + exterior[en + 1][EXT_WC] == exterior[st][EXT]) {
            computed.base_ctds[st] = CTD_FLUSH_COAX_WITH_NEXT;
            computed.base_ctds[en + 1] = CTD_FLUSH_COAX_WITH_PREV;
            q.emplace(st, en, DP_P);
            q.emplace(en + 1, -1, EXT_WC);
            goto loopend;
          }
          if ((enrb == G || enrb == U) && base00 + em.stack[enb][enrb][enrb ^ 1][stb] +
              exterior[en + 1][EXT_GU] == exterior[st][EXT]) {
            computed.base_ctds[st] = CTD_FLUSH_COAX_WITH_NEXT;
            computed.base_ctds[en + 1] = CTD_FLUSH_COAX_WITH_PREV;
            q.emplace(st, en, DP_P);
            q.emplace(en + 1, -1, EXT_GU);
            goto loopend;
          }
        }
      }
    } else {
      const auto stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];
      if (a == DP_P) {
        // It's paired, so add it to the folding.
        computed.s.p[st] = en;
        computed.s.p[en] = st;

        // Following largely matches the above DP so look up there for comments.
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (arr[ist][ien][DP_P] < CAP_E) {
              const auto val = em.TwoLoop(r, st, en, ist, ien) + arr[ist][ien][DP_P];
              if (val == arr[st][en][DP_P]) {
                q.emplace(ist, ien, DP_P);
                goto loopend;
              }
            }
          }
        }

        const auto base_branch_cost = em.AuGuPenalty(stb, enb) + em.multiloop_hack_a + em.multiloop_hack_b;
        // (<   ><    >)
        if (base_branch_cost + arr[st + 1][en - 1][DP_U2] == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_UNUSED;
          q.emplace(st + 1, en - 1, DP_U2);
          goto loopend;
        }
        // (3<   ><   >) 3'
        if (base_branch_cost + arr[st + 2][en - 1][DP_U2] + em.dangle3[stb][st1b][enb] == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_3_DANGLE;
          q.emplace(st + 2, en - 1, DP_U2);
          goto loopend;
        }
        // (<   ><   >5) 5'
        if (base_branch_cost + arr[st + 1][en - 2][DP_U2] + em.dangle5[stb][en1b][enb] == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_5_DANGLE;
          q.emplace(st + 1, en - 2, DP_U2);
          goto loopend;
        }
        // (.<   ><   >.) Terminal mismatch
        if (base_branch_cost + arr[st + 2][en - 2][DP_U2] + em.terminal[stb][st1b][en1b][enb] == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_TERMINAL_MISMATCH;
          q.emplace(st + 2, en - 2, DP_U2);
          goto loopend;
        }

        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          const base_t pl1b = r[piv - 1], plb = r[piv], prb = r[piv + 1], pr1b = r[piv + 2];

          // (.(   )   .) Left outer coax - P
          const auto outer_coax = em.MismatchCoaxial(stb, st1b, en1b, enb);
          if (base_branch_cost + arr[st + 2][piv][DP_P] + em.multiloop_hack_b +
              em.AuGuPenalty(st2b, plb) + arr[piv + 1][en - 2][DP_U] + outer_coax == arr[st][en][DP_P]) {
            computed.base_ctds[en] = CTD_LEFT_MISMATCH_COAX_WITH_NEXT;
            computed.base_ctds[st + 2] = CTD_LEFT_MISMATCH_COAX_WITH_PREV;
            q.emplace(st + 2, piv, DP_P);
            q.emplace(piv + 1, en - 2, DP_U);
            goto loopend;
          }
          // (.   (   ).) Right outer coax
          if (base_branch_cost + arr[st + 2][piv][DP_U] + em.multiloop_hack_b +
              em.AuGuPenalty(prb, en2b) + arr[piv + 1][en - 2][DP_P] + outer_coax == arr[st][en][DP_P]) {
            computed.base_ctds[en] = CTD_RIGHT_MISMATCH_COAX_WITH_PREV;
            computed.base_ctds[piv + 1] = CTD_RIGHT_MISMATCH_COAX_WITH_NEXT;
            q.emplace(st + 2, piv, DP_U);
            q.emplace(piv + 1, en - 2, DP_P);
            goto loopend;
          }

          // (.(   ).   ) Left right coax
          if (base_branch_cost + arr[st + 2][piv - 1][DP_P] + em.multiloop_hack_b +
              em.AuGuPenalty(st2b, pl1b) + arr[piv + 1][en - 1][DP_U] +
              em.MismatchCoaxial(pl1b, plb, st1b, st2b) == arr[st][en][DP_P]) {
            computed.base_ctds[en] = CTD_RIGHT_MISMATCH_COAX_WITH_NEXT;
            computed.base_ctds[st + 2] = CTD_RIGHT_MISMATCH_COAX_WITH_PREV;
            q.emplace(st + 2, piv - 1, DP_P);
            q.emplace(piv + 1, en - 1, DP_U);
            goto loopend;
          }
          // (   .(   ).) Right left coax
          if (base_branch_cost + arr[st + 1][piv][DP_U] + em.multiloop_hack_b +
              em.AuGuPenalty(pr1b, en2b) + arr[piv + 2][en - 2][DP_P] +
              em.MismatchCoaxial(en2b, en1b, prb, pr1b) == arr[st][en][DP_P]) {
            computed.base_ctds[en] = CTD_LEFT_MISMATCH_COAX_WITH_PREV;
            computed.base_ctds[piv + 2] = CTD_LEFT_MISMATCH_COAX_WITH_NEXT;
            q.emplace(st + 1, piv, DP_U);
            q.emplace(piv + 2, en - 2, DP_P);
            goto loopend;
          }

          // ((   )   ) Left flush coax
          if (base_branch_cost + arr[st + 1][piv][DP_P] +
              em.multiloop_hack_b + em.AuGuPenalty(st1b, plb) +
              arr[piv + 1][en - 1][DP_U] + em.stack[stb][st1b][plb][enb] == arr[st][en][DP_P]) {
            computed.base_ctds[en] = CTD_FLUSH_COAX_WITH_NEXT;
            computed.base_ctds[st + 1] = CTD_FLUSH_COAX_WITH_PREV;
            q.emplace(st + 1, piv, DP_P);
            q.emplace(piv + 1, en - 1, DP_U);
            goto loopend;
          }
          // (   (   )) Right flush coax
          if (base_branch_cost + arr[st + 1][piv][DP_U] +
              em.multiloop_hack_b + em.AuGuPenalty(prb, en1b) +
              arr[piv + 1][en - 1][DP_P] + em.stack[stb][prb][en1b][enb] == arr[st][en][DP_P]) {
            computed.base_ctds[en] = CTD_FLUSH_COAX_WITH_PREV;
            computed.base_ctds[piv + 1] = CTD_FLUSH_COAX_WITH_NEXT;
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
      if (st + 1 < en && (a == DP_U || a == DP_U2) && arr[st + 1][en][a] == arr[st][en][a]) {
        q.emplace(st + 1, en, a);
        goto loopend;
      }

      // Pair here.
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = r[piv], pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const auto base00 = arr[st][piv][DP_P] + em.AuGuPenalty(stb, pb) + em.multiloop_hack_b;
        const auto base01 = arr[st][piv - 1][DP_P] + em.AuGuPenalty(stb, pl1b) + em.multiloop_hack_b;
        const auto base10 = arr[st + 1][piv][DP_P] + em.AuGuPenalty(st1b, pb) + em.multiloop_hack_b;
        const auto base11 = arr[st + 1][piv - 1][DP_P] + em.AuGuPenalty(st1b, pl1b) + em.multiloop_hack_b;

        // Min is for either placing another unpaired or leaving it as nothing.
        // If we're at U2, don't allow leaving as nothing.
        auto right_unpaired = arr[piv + 1][en][DP_U];
        if (a != DP_U2)
          right_unpaired = std::min(right_unpaired, 0);

        // Check a == U_RCOAX:
        // (   ).<( ** ). > Right coax backward
        if (a == DP_U_RCOAX) {
          if (st > 0 && base01 + em.MismatchCoaxial(
              pl1b, pb, r[st - 1], stb) + right_unpaired == arr[st][en][DP_U_RCOAX]) {
            // Ctds were already set from the recurrence that called this.
            q.emplace(st, piv - 1, DP_P);
            if (right_unpaired)
              q.emplace(piv + 1, en, DP_U);
            goto loopend;
          }
          continue;
        }

        // (   )<   > - U, U2, U_WC?, U_GU?
        if (base00 + right_unpaired == arr[st][en][a] &&
            (a != DP_U_WC || IsWatsonCrick(stb, pb)) &&
            (a != DP_U_GU || IsGu(stb, pb))) {
          // If U_WC, or U_GU, we were involved in some sort of coaxial stack previously, and were already set.
          if (a != DP_U_WC && a != DP_U_GU)
            computed.base_ctds[st] = CTD_UNUSED;
          q.emplace(st, piv, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }

        // The rest of the cases are for U and U2.
        if (a != DP_U && a != DP_U2)
          continue;

        // (   )3<   > 3' - U, U2
        if (base01 + em.dangle3[pl1b][pb][stb] + right_unpaired == arr[st][en][a]) {
          computed.base_ctds[st] = CTD_3_DANGLE;
          q.emplace(st, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // 5(   )<   > 5' - U, U2
        if (base10 + em.dangle5[pb][stb][st1b] + right_unpaired == arr[st][en][a]) {
          computed.base_ctds[st + 1] = CTD_5_DANGLE;
          q.emplace(st + 1, piv, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch - U, U2
        if (base11 + em.terminal[pl1b][pb][stb][st1b] + right_unpaired == arr[st][en][a]) {
          computed.base_ctds[st + 1] = CTD_TERMINAL_MISMATCH;
          q.emplace(st + 1, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<(   ) > Left coax - U, U2
        auto val = base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b);
        if (val + arr[piv + 1][en][DP_U_WC] == arr[st][en][a]) {
          computed.base_ctds[st + 1] = CTD_LEFT_MISMATCH_COAX_WITH_NEXT;
          computed.base_ctds[piv + 1] = CTD_LEFT_MISMATCH_COAX_WITH_PREV;
          q.emplace(st + 1, piv - 1, DP_P);
          q.emplace(piv + 1, en, DP_U_WC);
          goto loopend;
        }
        if (val + arr[piv + 1][en][DP_U_GU] == arr[st][en][a]) {
          computed.base_ctds[st + 1] = CTD_LEFT_MISMATCH_COAX_WITH_NEXT;
          computed.base_ctds[piv + 1] = CTD_LEFT_MISMATCH_COAX_WITH_PREV;
          q.emplace(st + 1, piv - 1, DP_P);
          q.emplace(piv + 1, en, DP_U_GU);
          goto loopend;
        }

        // (   ).<(   ). > Right coax forward - U, U2
        if (base01 + arr[piv + 1][en][DP_U_RCOAX] == arr[st][en][a]) {
          computed.base_ctds[st] = CTD_RIGHT_MISMATCH_COAX_WITH_NEXT;
          computed.base_ctds[piv + 1] = CTD_RIGHT_MISMATCH_COAX_WITH_PREV;
          q.emplace(st, piv - 1, DP_P);
          q.emplace(piv + 1, en, DP_U_RCOAX);
          goto loopend;
        }

        // There has to be remaining bases to even have a chance at these cases.
        if (piv < en) {
          const auto pr1b = r[piv + 1];
          // (   )<(   ) > Flush coax - U, U2
          if (base00 + em.stack[pb][pr1b][pr1b ^ 3][stb] + arr[piv + 1][en][DP_U_WC] == arr[st][en][a]) {
            computed.base_ctds[st] = CTD_FLUSH_COAX_WITH_NEXT;
            computed.base_ctds[piv + 1] = CTD_FLUSH_COAX_WITH_PREV;
            q.emplace(st, piv, DP_P);
            q.emplace(piv + 1, en, DP_U_WC);
            goto loopend;
          }
          if ((pr1b == G || pr1b == U) &&
              base00 + em.stack[pb][pr1b][pr1b ^ 1][stb] + arr[piv + 1][en][DP_U_GU] == arr[st][en][a]) {
            computed.base_ctds[st] = CTD_FLUSH_COAX_WITH_NEXT;
            computed.base_ctds[piv + 1] = CTD_FLUSH_COAX_WITH_PREV;
            q.emplace(st, piv, DP_P);
            q.emplace(piv + 1, en, DP_U_GU);
            goto loopend;
          }
        }
      }
    }

    loopend:;
  }
  return computed;
}


}
}
