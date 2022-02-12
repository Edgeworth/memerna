// Copyright 2016 Eliot Courtney.
#include <algorithm>

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "compute/mfe/mfe.h"
#include "model/base.h"
#include "model/model.h"
#include "model/primary.h"
#include "util/array.h"

namespace mrna::mfe {

#define UPDATE_CACHE(a, value)                                          \
  do {                                                                  \
    Energy macro_upd_value_ = (value);                                  \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < dp[st][en][a]) { \
      dp[st][en][a] = macro_upd_value_;                                 \
    }                                                                   \
  } while (0)

DpArray ComputeTables0(const Primary& r, const energy::EnergyModelPtr& em) {
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  const int N = static_cast<int>(r.size());
  auto dp = DpArray(r.size() + 1, MAX_E);

  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const Base stb = r[st];
      const Base st1b = r[st + 1];
      const Base st2b = r[st + 2];
      const Base enb = r[en];
      const Base en1b = r[en - 1];
      const Base en2b = r[en - 2];

      // TODO: check lonely pairs
      if (em->CanPair(r, st, en)) {
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (dp[ist][ien][DP_P] < CAP_E)
              UPDATE_CACHE(DP_P, em->TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        UPDATE_CACHE(DP_P, em->Hairpin(r, st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost =
            em->AuGuPenalty(stb, enb) + em->multiloop_hack_a + em->multiloop_hack_b;

        // (<   ><   >)
        UPDATE_CACHE(DP_P, base_branch_cost + dp[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        UPDATE_CACHE(
            DP_P, base_branch_cost + dp[st + 2][en - 1][DP_U2] + em->dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        UPDATE_CACHE(
            DP_P, base_branch_cost + dp[st + 1][en - 2][DP_U2] + em->dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        UPDATE_CACHE(DP_P,
            base_branch_cost + dp[st + 2][en - 2][DP_U2] + em->terminal[stb][st1b][en1b][enb]);

        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          Base pl1b = r[piv - 1];
          Base plb = r[piv];
          Base prb = r[piv + 1];
          Base pr1b = r[piv + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          const auto outer_coax = em->MismatchCoaxial(stb, st1b, en1b, enb);
          UPDATE_CACHE(DP_P,
              base_branch_cost + dp[st + 2][piv][DP_P] + em->multiloop_hack_b +
                  em->AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax);
          // (.   (   ).) Right outer coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + dp[st + 2][piv][DP_U] + em->multiloop_hack_b +
                  em->AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax);

          // (.(   ).   ) Left right coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + dp[st + 2][piv - 1][DP_P] + em->multiloop_hack_b +
                  em->AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
                  em->MismatchCoaxial(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + dp[st + 1][piv][DP_U] + em->multiloop_hack_b +
                  em->AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
                  em->MismatchCoaxial(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + dp[st + 1][piv][DP_P] + em->multiloop_hack_b +
                  em->AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] +
                  em->stack[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + dp[st + 1][piv][DP_U] + em->multiloop_hack_b +
                  em->AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] +
                  em->stack[stb][prb][en1b][enb]);
        }
      }

      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        UPDATE_CACHE(DP_U, dp[st + 1][en][DP_U]);
        UPDATE_CACHE(DP_U2, dp[st + 1][en][DP_U2]);
      }
      // Pair here.
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = r[piv];
        const auto pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const auto base00 = dp[st][piv][DP_P] + em->AuGuPenalty(stb, pb) + em->multiloop_hack_b;
        const auto base01 =
            dp[st][piv - 1][DP_P] + em->AuGuPenalty(stb, pl1b) + em->multiloop_hack_b;
        const auto base10 =
            dp[st + 1][piv][DP_P] + em->AuGuPenalty(st1b, pb) + em->multiloop_hack_b;
        const auto base11 =
            dp[st + 1][piv - 1][DP_P] + em->AuGuPenalty(st1b, pl1b) + em->multiloop_hack_b;
        // Min is for either placing another unpaired or leaving it as nothing.
        const auto right_unpaired = std::min(dp[piv + 1][en][DP_U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        UPDATE_CACHE(DP_U2, base00 + dp[piv + 1][en][DP_U]);
        auto val = base00 + right_unpaired;
        UPDATE_CACHE(DP_U, val);
        if (IsGuPair(stb, pb))
          UPDATE_CACHE(DP_U_GU, val);
        else
          UPDATE_CACHE(DP_U_WC, val);

        // (   )3<   > 3' - U
        UPDATE_CACHE(DP_U, base01 + em->dangle3[pl1b][pb][stb] + right_unpaired);
        UPDATE_CACHE(DP_U2, base01 + em->dangle3[pl1b][pb][stb] + dp[piv + 1][en][DP_U]);
        // 5(   )<   > 5' - U
        UPDATE_CACHE(DP_U, base10 + em->dangle5[pb][stb][st1b] + right_unpaired);
        UPDATE_CACHE(DP_U2, base10 + em->dangle5[pb][stb][st1b] + dp[piv + 1][en][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        UPDATE_CACHE(DP_U, base11 + em->terminal[pl1b][pb][stb][st1b] + right_unpaired);
        UPDATE_CACHE(DP_U2, base11 + em->terminal[pl1b][pb][stb][st1b] + dp[piv + 1][en][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + em->MismatchCoaxial(pl1b, pb, stb, st1b) +
            std::min(dp[piv + 1][en][DP_U_WC], dp[piv + 1][en][DP_U_GU]);
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);

        // (   )<.(   ). > Right coax forward and backward
        val = base00 + dp[piv + 1][en][DP_U_RCOAX];
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);
        UPDATE_CACHE(
            DP_U_RCOAX, base11 + em->MismatchCoaxial(pl1b, pb, stb, st1b) + right_unpaired);

        // (   )(<   ) > Flush coax - U
        val = base01 + em->stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC];
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);
        if (IsGu(pb)) {
          val = base01 + em->stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU];
          UPDATE_CACHE(DP_U, val);
          UPDATE_CACHE(DP_U2, val);
        }
      }
    }
  }

  return dp;
}

#undef UPDATE_CACHE

}  // namespace mrna::mfe
