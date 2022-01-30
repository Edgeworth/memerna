// Copyright 2016 Eliot Courtney.
#include <algorithm>

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "compute/energy/precomp.h"
#include "compute/mfe/mfe.h"
#include "model/base.h"
#include "model/model.h"
#include "model/primary.h"
#include "util/array.h"

namespace mrna::mfe {

DpArray ComputeTables1(const Primary& r, const energy::EnergyModel& em) {
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  const int N = static_cast<int>(r.size());
  const energy::Precomp pc(Primary(r), em);
  auto dp = DpArray(r.size() + 1, MAX_E);

  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const Base stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1],
                 en2b = r[en - 2];

      // TODO: check lonely pairs
      if (em.CanPair(r, st, en)) {
        Energy p_min = MAX_E;
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (dp[ist][ien][DP_P] < CAP_E)
              p_min = std::min(p_min, pc.TwoLoop(st, en, ist, ien) + dp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        p_min = std::min(p_min, em.Hairpin(r, st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost = pc.augubranch[stb][enb] + em.multiloop_hack_a;

        // (<   ><   >)
        p_min = std::min(p_min, base_branch_cost + dp[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        p_min = std::min(
            p_min, base_branch_cost + dp[st + 2][en - 1][DP_U2] + em.dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        p_min = std::min(
            p_min, base_branch_cost + dp[st + 1][en - 2][DP_U2] + em.dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        p_min = std::min(p_min,
            base_branch_cost + dp[st + 2][en - 2][DP_U2] + em.terminal[stb][st1b][en1b][enb]);

        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          Base pl1b = r[piv - 1], plb = r[piv], prb = r[piv + 1], pr1b = r[piv + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          const auto outer_coax = em.MismatchCoaxial(stb, st1b, en1b, enb);
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 2][piv][DP_P] + pc.augubranch[st2b][plb] +
                  dp[piv + 1][en - 2][DP_U] + outer_coax);
          // (.   (   ).) Right outer coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 2][piv][DP_U] + pc.augubranch[prb][en2b] +
                  dp[piv + 1][en - 2][DP_P] + outer_coax);

          // (.(   ).   ) Left right coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 2][piv - 1][DP_P] + pc.augubranch[st2b][pl1b] +
                  dp[piv + 1][en - 1][DP_U] + em.MismatchCoaxial(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 1][piv][DP_U] + pc.augubranch[pr1b][en2b] +
                  dp[piv + 2][en - 2][DP_P] + em.MismatchCoaxial(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 1][piv][DP_P] + pc.augubranch[st1b][plb] +
                  dp[piv + 1][en - 1][DP_U] + em.stack[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 1][piv][DP_U] + pc.augubranch[prb][en1b] +
                  dp[piv + 1][en - 1][DP_P] + em.stack[stb][prb][en1b][enb]);
        }

        dp[st][en][DP_P] = p_min;
      }
      Energy u_min = MAX_E, u2_min = MAX_E, rcoax_min = MAX_E, wc_min = MAX_E, gu_min = MAX_E;
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u_min = std::min(u_min, dp[st + 1][en][DP_U]);
        u2_min = std::min(u2_min, dp[st + 1][en][DP_U2]);
      }
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = r[piv], pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const auto base00 = dp[st][piv][DP_P] + pc.augubranch[stb][pb];
        const auto base01 = dp[st][piv - 1][DP_P] + pc.augubranch[stb][pl1b];
        const auto base10 = dp[st + 1][piv][DP_P] + pc.augubranch[st1b][pb];
        const auto base11 = dp[st + 1][piv - 1][DP_P] + pc.augubranch[st1b][pl1b];
        // Min is for either placing another unpaired or leaving it as nothing.
        const auto right_unpaired = std::min(dp[piv + 1][en][DP_U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        u2_min = std::min(u2_min, base00 + dp[piv + 1][en][DP_U]);
        auto val = base00 + right_unpaired;
        u_min = std::min(u_min, val);
        if (IsGu(stb, pb))
          gu_min = std::min(gu_min, val);
        else
          wc_min = std::min(wc_min, val);

        // (   )3<   > 3' - U
        u_min = std::min(u_min, base01 + em.dangle3[pl1b][pb][stb] + right_unpaired);
        u2_min = std::min(u2_min, base01 + em.dangle3[pl1b][pb][stb] + dp[piv + 1][en][DP_U]);
        // 5(   )<   > 5' - U
        u_min = std::min(u_min, base10 + em.dangle5[pb][stb][st1b] + right_unpaired);
        u2_min = std::min(u2_min, base10 + em.dangle5[pb][stb][st1b] + dp[piv + 1][en][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        u_min = std::min(u_min, base11 + em.terminal[pl1b][pb][stb][st1b] + right_unpaired);
        u2_min =
            std::min(u2_min, base11 + em.terminal[pl1b][pb][stb][st1b] + dp[piv + 1][en][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b) +
            std::min(dp[piv + 1][en][DP_U_WC], dp[piv + 1][en][DP_U_GU]);
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);

        // (   )<.(   ). > Right coax forward and backward
        val = base00 + dp[piv + 1][en][DP_U_RCOAX];
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
        rcoax_min =
            std::min(rcoax_min, base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b) + right_unpaired);

        // (   )(<   ) > Flush coax - U
        val = base01 + em.stack[pl1b][pb][pb ^ 3][stb] + dp[piv][en][DP_U_WC];
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
        if (pb == G || pb == U) {
          val = base01 + em.stack[pl1b][pb][pb ^ 1][stb] + dp[piv][en][DP_U_GU];
          u_min = std::min(u_min, val);
          u2_min = std::min(u2_min, val);
        }
      }

      dp[st][en][DP_U] = u_min;
      dp[st][en][DP_U2] = u2_min;
      dp[st][en][DP_U_WC] = wc_min;
      dp[st][en][DP_U_GU] = gu_min;
      dp[st][en][DP_U_RCOAX] = rcoax_min;
    }
  }

  return dp;
}

}  // namespace mrna::mfe
