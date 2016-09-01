#include "traceback.h"
#include "fold/fold_internal.h"

namespace memerna {
namespace fold {

using namespace constants;
using namespace energy;
using namespace internal;

#define UPDATE_EXT(a, na, pst, pen, value) \
  do { \
    energy_t macro_upd_value_ = (value) + exterior[en + 1][na]; \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < exterior[st][a]) { \
      exterior[st][a] = macro_upd_value_; \
      exterior_sts[st][a] = std::make_tuple(pst, pen, en + 1, na); \
    } \
  } while (0)

array2d_t<energy_t, EXT_SIZE> TraceExterior(const primary_t& r,
    const array3d_t<energy_t, DP_SIZE>& arr, traceback_stack_t& q) {
  int N = int(r.size());
  // Exterior loop calculation. There can be no paired base on exterior[en].
  array2d_t<energy_t, EXT_SIZE> exterior(std::size_t(N + 1));
  exterior[N][EXT] = 0;
  // Holds: st, en, nst, na - nst and na for next index into itself; st, en for the paired.
  array2d_t<std::tuple<int, int, int, int>, EXT_SIZE> exterior_sts(std::size_t(N), 0xFF);
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    exterior[st][EXT] = exterior[st + 1][EXT];
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      auto stb = r[st], st1b = r[st + 1], enb = r[en], en1b = r[en - 1];

      auto base00 = arr[st][en][DP_P] + AuGuPenalty(stb, enb);
      auto base01 = arr[st][en - 1][DP_P] + AuGuPenalty(stb, en1b);
      auto base10 = arr[st + 1][en][DP_P] + AuGuPenalty(st1b, enb);
      auto base11 = arr[st + 1][en - 1][DP_P] + AuGuPenalty(st1b, en1b);

      // (   )<   >
      UPDATE_EXT(EXT, EXT, st, en, base00);
      if (IsGu(stb, enb))
        UPDATE_EXT(EXT_GU, EXT, st, en, base00);
      else
        UPDATE_EXT(EXT_WC, EXT, st, en, base00);

      // (   )3<   > 3'
      UPDATE_EXT(EXT, EXT, st, en - 1, base01 + g_dangle3[en1b][enb][stb]);
      // 5(   )<   > 5'
      UPDATE_EXT(EXT, EXT, st + 1, en, base10 + g_dangle5[enb][stb][st1b]);
      // .(   ).<   > Terminal mismatch
      UPDATE_EXT(EXT, EXT, st + 1, en - 1, base11 + g_terminal[en1b][enb][stb][st1b]);
      // .(   ).<(   ) > Left coax  x
      auto val = base11 + MismatchCoaxial(en1b, enb, stb, st1b);
      UPDATE_EXT(EXT, EXT_GU, st + 1, en - 1, val);
      UPDATE_EXT(EXT, EXT_WC, st + 1, en - 1, val);

      // (   ).<(   ). > Right coax forward
      UPDATE_EXT(EXT, EXT_RCOAX, st, en - 1, base01);
      // (   ).<( * ). > Right coax backward
      if (st > 0)
        UPDATE_EXT(EXT_RCOAX, EXT, st, en - 1, base01 + MismatchCoaxial(
            en1b, enb, r[st - 1], stb));

      if (en < N - 1) {
        // (   )<(   ) > Flush coax
        auto enrb = r[en + 1];
        UPDATE_EXT(EXT, EXT_WC, st, en, base00 + g_stack[enb][enrb][enrb ^ 3][stb]);
        if (enrb == G || enrb == U)
          UPDATE_EXT(EXT, EXT_GU, st, en, base00 + g_stack[enb][enrb][enrb ^ 1][stb]);
      }
    }
  }

  // Backtrace.
  // sz, st, paired
  int ext_st = 0, ext_a = EXT;
  while (ext_st < N) {
    int psz, pst, nst, na;
    std::tie(psz, pst, nst, na) = exterior_sts[ext_st][ext_a];
    if (nst == -1) {
      ++ext_st;
      continue;
    }
    q.emplace(psz, pst, DP_P);
    ext_st = nst;
    ext_a = na;
  }
  return exterior;
}

#undef UPDATE_EXT

computed_t TraceStructure(const primary_t& r, const array3d_t<energy_t, DP_SIZE>& arr,
    const array2d_t<energy_t, EXT_SIZE>& exterior, traceback_stack_t& q) {
  computed_t computed(r);
  computed.energy = exterior[0][EXT];
  while (!q.empty()) {
    int st, en, a;
    std::tie(st, en, a) = q.top();
    //printf("%d %d %d\n", st, en, a);
    auto stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];
    q.pop();
    if (a == DP_P) {
      // It's paired, so add it to the folding.
      computed.s.p[st] = en;
      computed.s.p[en] = st;

      // Following largely matches the above DP so look up there for comments.
      int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
      for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
        for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
          if (arr[ist][ien][DP_P] < CAP_E) {
            auto val = TwoLoop(r, st, en, ist, ien) + arr[ist][ien][DP_P];
            if (val == arr[st][en][DP_P]) {
              q.emplace(ist, ien, DP_P);
              goto loopend;
            }
          }
        }
      }

      auto base_branch_cost = AuGuPenalty(stb, enb) + g_multiloop_hack_a + g_multiloop_hack_b;

      // (<   ><    >)
      if (base_branch_cost + arr[st + 1][en - 1][DP_U2] == arr[st][en][DP_P]) {
        computed.base_ctds[en] = CTD_UNUSED;
        q.emplace(st + 1, en - 1, DP_U2);
        goto loopend;
      }
      // (3<   ><   >) 3'
      if (base_branch_cost + arr[st + 2][en - 1][DP_U2] + g_dangle3[stb][st1b][enb] == arr[st][en][DP_P]) {
        computed.base_ctds[en] = CTD_3_DANGLE;
        q.emplace(st + 2, en - 1, DP_U2);
        goto loopend;
      }
      // (<   ><   >5) 5'
      if (base_branch_cost + arr[st + 1][en - 2][DP_U2] + g_dangle5[stb][en1b][enb] == arr[st][en][DP_P]) {
        computed.base_ctds[en] = CTD_5_DANGLE;
        q.emplace(st + 1, en - 2, DP_U2);
        goto loopend;
      }
      // (.<   ><   >.) Terminal mismatch
      if (base_branch_cost + arr[st + 2][en - 2][DP_U2] + g_terminal[stb][st1b][en1b][enb] == arr[st][en][DP_P]) {
        computed.base_ctds[en] = CTD_TERMINAL_MISMATCH;
        q.emplace(st + 2, en - 2, DP_U2);
        goto loopend;
      }

      for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
        base_t pl1b = r[piv - 1], plb = r[piv], prb = r[piv + 1], pr1b = r[piv + 2];

        // (.(   )   .) Left outer coax - P
        auto outer_coax = MismatchCoaxial(stb, st1b, en1b, enb);
        if (base_branch_cost + arr[st + 2][piv][DP_P] + g_multiloop_hack_b +
            AuGuPenalty(st2b, plb) + arr[piv + 1][en - 2][DP_U] + outer_coax == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_LEFT_MISMATCH_COAX_WITH_NEXT;
          computed.base_ctds[st + 2] = CTD_LEFT_MISMATCH_COAX_WITH_PREV;
          q.emplace(st + 2, piv, DP_P);
          q.emplace(piv + 1, en - 2, DP_U);
          goto loopend;
        }
        // (.   (   ).) Right outer coax
        if (base_branch_cost + arr[st + 2][piv][DP_U] + g_multiloop_hack_b +
            AuGuPenalty(prb, en2b) + arr[piv + 1][en - 2][DP_P] + outer_coax == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_RIGHT_MISMATCH_COAX_WITH_PREV;
          computed.base_ctds[piv + 1] = CTD_RIGHT_MISMATCH_COAX_WITH_NEXT;
          q.emplace(st + 2, piv, DP_U);
          q.emplace(piv + 1, en - 2, DP_P);
          goto loopend;
        }

        // (.(   ).   ) Left right coax
        if (base_branch_cost + arr[st + 2][piv - 1][DP_P] + g_multiloop_hack_b +
            AuGuPenalty(st2b, pl1b) + arr[piv + 1][en - 1][DP_U] +
            MismatchCoaxial(pl1b, plb, st1b, st2b) == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_RIGHT_MISMATCH_COAX_WITH_NEXT;
          computed.base_ctds[st + 2] = CTD_RIGHT_MISMATCH_COAX_WITH_PREV;
          q.emplace(st + 2, piv - 1, DP_P);
          q.emplace(piv + 1, en - 1, DP_U);
          goto loopend;
        }
        // (   .(   ).) Right left coax
        if (base_branch_cost + arr[st + 1][piv][DP_U] + g_multiloop_hack_b +
            AuGuPenalty(pr1b, en2b) + arr[piv + 2][en - 2][DP_P] +
            MismatchCoaxial(en2b, en1b, prb, pr1b) == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_LEFT_MISMATCH_COAX_WITH_PREV;
          computed.base_ctds[piv + 2] = CTD_LEFT_MISMATCH_COAX_WITH_NEXT;
          q.emplace(st + 1, piv, DP_U);
          q.emplace(piv + 2, en - 2, DP_P);
          goto loopend;
        }

        // ((   )   ) Left flush coax
        if (base_branch_cost + arr[st + 1][piv][DP_P] +
            g_multiloop_hack_b + AuGuPenalty(st1b, plb) +
            arr[piv + 1][en - 1][DP_U] + g_stack[stb][st1b][plb][enb] == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_FLUSH_COAX_WITH_NEXT;
          computed.base_ctds[st + 1] = CTD_FLUSH_COAX_WITH_PREV;
          q.emplace(st + 1, piv, DP_P);
          q.emplace(piv + 1, en - 1, DP_U);
          goto loopend;
        }
        // (   (   )) Right flush coax
        if (base_branch_cost + arr[st + 1][piv][DP_U] +
            g_multiloop_hack_b + AuGuPenalty(prb, en1b) +
            arr[piv + 1][en - 1][DP_P] + g_stack[stb][prb][en1b][enb] == arr[st][en][DP_P]) {
          computed.base_ctds[en] = CTD_FLUSH_COAX_WITH_PREV;
          computed.base_ctds[piv + 1] = CTD_FLUSH_COAX_WITH_NEXT;
          q.emplace(st + 1, piv, DP_U);
          q.emplace(piv + 1, en - 1, DP_P);
          goto loopend;
        }
      }
    } else {
      // Left unpaired. Either DP_U or DP_U2.
      if (st + 1 < en && (a == DP_U || a == DP_U2) && arr[st + 1][en][a] == arr[st][en][a]) {
        q.emplace(st + 1, en, a);
        goto loopend;
      }

      // Pair here.
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        auto pb = r[piv], pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[st][piv][DP_P] + AuGuPenalty(stb, pb) + g_multiloop_hack_b;
        auto base01 = arr[st][piv - 1][DP_P] + AuGuPenalty(stb, pl1b) + g_multiloop_hack_b;
        auto base10 = arr[st + 1][piv][DP_P] + AuGuPenalty(st1b, pb) + g_multiloop_hack_b;
        auto base11 = arr[st + 1][piv - 1][DP_P] + AuGuPenalty(st1b, pl1b) + g_multiloop_hack_b;

        // Min is for either placing another unpaired or leaving it as nothing.
        // If we're at U2, don't allow leaving as nothing.
        auto right_unpaired = arr[piv + 1][en][DP_U];
        if (a != DP_U2)
          right_unpaired = std::min(right_unpaired, 0);

        // Check a == U_RCOAX:
        // (   ).<( ** ). > Right coax backward
        if (a == DP_U_RCOAX) {
          if (st > 0 && base01 + MismatchCoaxial(
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
        if (base01 + g_dangle3[pl1b][pb][stb] + right_unpaired == arr[st][en][a]) {
          computed.base_ctds[st] = CTD_3_DANGLE;
          q.emplace(st, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // 5(   )<   > 5' - U, U2
        if (base10 + g_dangle5[pb][stb][st1b] + right_unpaired == arr[st][en][a]) {
          computed.base_ctds[st + 1] = CTD_5_DANGLE;
          q.emplace(st + 1, piv, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch - U, U2
        if (base11 + g_terminal[pl1b][pb][stb][st1b] + right_unpaired == arr[st][en][a]) {
          computed.base_ctds[st + 1] = CTD_TERMINAL_MISMATCH;
          q.emplace(st + 1, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<(   ) > Left coax - U, U2
        auto val = base11 + MismatchCoaxial(pl1b, pb, stb, st1b);
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
          auto pr1b = r[piv + 1];
          // (   )<(   ) > Flush coax - U, U2
          if (base00 + g_stack[pb][pr1b][pr1b ^ 3][stb] + arr[piv + 1][en][DP_U_WC] == arr[st][en][a]) {
            computed.base_ctds[st] = CTD_FLUSH_COAX_WITH_NEXT;
            computed.base_ctds[piv + 1] = CTD_FLUSH_COAX_WITH_PREV;
            q.emplace(st, piv, DP_P);
            q.emplace(piv + 1, en, DP_U_WC);
            goto loopend;
          }
          if ((pr1b == G || pr1b == U) &&
              base00 + g_stack[pb][pr1b][pr1b ^ 1][stb] + arr[piv + 1][en][DP_U_GU] == arr[st][en][a]) {
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
