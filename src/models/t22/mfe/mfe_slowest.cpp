// Copyright 2023 E.
#include <algorithm>
#include <compare>
#include <memory>

#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "models/t04/mfe/dp.h"
#include "models/t22/energy/model.h"
#include "models/t22/mfe/mfe.h"

namespace mrna::md::t22 {

using t04::DP_P;
using t04::DP_U;
using t04::DP_U2;
using t04::DP_U_GU;
using t04::DP_U_RC;
using t04::DP_U_WC;
using t04::DpArray;

// TODO(0): Implement. Any way to generalise/pull out code?
void MfeSlowest(const Primary& r, const Model::Ptr& em, DpState& state) {
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  const int N = static_cast<int>(r.size());
  state.t04.dp = DpArray(r.size() + 1, MAX_E);
  state.penult = Array3D(r.size() + 1, MAX_E);
  auto& dp = state.t04.dp;
  auto& penult = state.penult;

  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const Base stb = r[st];
      const Base st1b = r[st + 1];
      const Base st2b = r[st + 2];
      const Base enb = r[en];
      const Base en1b = r[en - 1];
      const Base en2b = r[en - 2];

      if (em->CanPair(r, st, en)) {
        Energy p_min = MAX_E;

        {
          const int max_stack = en - st - HAIRPIN_MIN_SZ + 1;
          // Try all stacks of each length, with or without a 1 nuc bulge loop intercedeing.
          const Energy bulge_left = em->Bulge(r, st, en, st + 2, en - 1);
          const Energy bulge_right = em->Bulge(r, st, en, st + 1, en - 2);
          const Energy augu_penalty = em->AuGuPenalty(stb, enb);

          // Base case of length 1 is not exactly possible, but to simplify the
          // implementation it contains the penultimate penalty assuming this is
          // the ending inner pair, if (st, en) can pair.
          if (st > 0 && en < N - 1) {
            // Can either continue by adding a multiloop, adding a internal
            // loop, or ending with a hairpin.
            auto best = std::min(dp[st + 1][en - 1][DP_U2] + augu_penalty, em->Hairpin(r, st, en));

            const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
            for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
              for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
                // Try all internal loops. We don't check stacks or 1 nuc bulge loops.
                if (dp[ist][ien][DP_P] < CAP_E && ist - st + en - ien > 3)
                  best = std::min(best, em->TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P]);
              }
            }

            penult[st][en][1] = std::min(
                penult[st][en][1], em->penultimate_stack[r[st - 1]][stb][enb][r[en + 1]] + best);
          }

          // Try stems with a specific length.
          for (int length = 2; 2 * length <= max_stack; ++length) {
            // Update our DP at (st, en) - no bulge:
            auto none =
                penult[st + 1][en - 1][length - 1] + em->stack[r[st]][r[st + 1]][r[en - 1]][r[en]];
            penult[st][en][length] = std::min(penult[st][en][length], none);
            // Left bulge:
            auto left = penult[st + 2][en - 1][length - 1] + bulge_left;
            penult[st][en][length] = std::min(penult[st][en][length], left);
            // Right bulge:
            auto right = penult[st + 1][en - 2][length - 1] + bulge_right;
            penult[st][en][length] = std::min(penult[st][en][length], right);

            // Try each of the bulges for the paired array.
            p_min =
                std::min(p_min, none + em->penultimate_stack[en1b][enb][stb][st1b] + augu_penalty);
            p_min =
                std::min(p_min, left + em->penultimate_stack[en1b][enb][stb][st2b] + augu_penalty);
            p_min =
                std::min(p_min, right + em->penultimate_stack[en2b][enb][stb][st1b] + augu_penalty);
          }
        }

        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            // Try all internal loops. We don't check stacks or 1 nuc bulge loops.
            if (dp[ist][ien][DP_P] < CAP_E && ist - st + en - ien > 3)
              p_min = std::min(p_min, em->TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        p_min = std::min(p_min, em->Hairpin(r, st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost =
            em->AuGuPenalty(stb, enb) + em->multiloop_hack_a + em->multiloop_hack_b;

        // (<   ><   >)
        p_min = std::min(p_min, base_branch_cost + dp[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        p_min = std::min(
            p_min, base_branch_cost + dp[st + 2][en - 1][DP_U2] + em->dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        p_min = std::min(
            p_min, base_branch_cost + dp[st + 1][en - 2][DP_U2] + em->dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        p_min = std::min(p_min,
            base_branch_cost + dp[st + 2][en - 2][DP_U2] + em->terminal[stb][st1b][en1b][enb]);

        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          const Base pl1b = r[piv - 1];
          const Base plb = r[piv];
          const Base prb = r[piv + 1];
          const Base pr1b = r[piv + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          const auto outer_coax = em->MismatchCoaxial(stb, st1b, en1b, enb);
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 2][piv][DP_P] + em->multiloop_hack_b +
                  em->AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax);
          // (.   (   ).) Right outer coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 2][piv][DP_U] + em->multiloop_hack_b +
                  em->AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax);

          // (.(   ).   ) Left inner coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 2][piv - 1][DP_P] + em->multiloop_hack_b +
                  em->AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
                  em->MismatchCoaxial(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right inner coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 1][piv][DP_U] + em->multiloop_hack_b +
                  em->AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
                  em->MismatchCoaxial(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 1][piv][DP_P] + em->multiloop_hack_b +
                  em->AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] +
                  em->stack[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 1][piv][DP_U] + em->multiloop_hack_b +
                  em->AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] +
                  em->stack[stb][prb][en1b][enb]);
        }

        dp[st][en][DP_P] = p_min;
      }
      Energy u_min = MAX_E;
      Energy u2_min = MAX_E;
      Energy rcoax_min = MAX_E;
      Energy wc_min = MAX_E;
      Energy gu_min = MAX_E;
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u_min = std::min(u_min, dp[st + 1][en][DP_U]);
        u2_min = std::min(u2_min, dp[st + 1][en][DP_U2]);
      }
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = r[piv];
        const auto pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the
        // right.
        const auto base00 = dp[st][piv][DP_P] + em->multiloop_hack_b + em->AuGuPenalty(stb, pb);
        const auto base01 =
            dp[st][piv - 1][DP_P] + em->multiloop_hack_b + em->AuGuPenalty(stb, pl1b);
        const auto base10 =
            dp[st + 1][piv][DP_P] + em->multiloop_hack_b + em->AuGuPenalty(st1b, pb);
        const auto base11 =
            dp[st + 1][piv - 1][DP_P] + em->multiloop_hack_b + em->AuGuPenalty(st1b, pl1b);
        // Min is for either placing another unpaired or leaving it as nothing.
        const auto right_unpaired = std::min(dp[piv + 1][en][DP_U], ZERO_E);

        // (   )<   > - U, U_WC?, U_GU?
        u2_min = std::min(u2_min, base00 + dp[piv + 1][en][DP_U]);
        auto val = base00 + right_unpaired;
        u_min = std::min(u_min, val);
        if (IsGuPair(stb, pb))
          gu_min = std::min(gu_min, val);
        else
          wc_min = std::min(wc_min, val);

        // (   )3<   > 3' - U
        u_min = std::min(u_min, base01 + em->dangle3[pl1b][pb][stb] + right_unpaired);
        u2_min = std::min(u2_min, base01 + em->dangle3[pl1b][pb][stb] + dp[piv + 1][en][DP_U]);
        // 5(   )<   > 5' - U
        u_min = std::min(u_min, base10 + em->dangle5[pb][stb][st1b] + right_unpaired);
        u2_min = std::min(u2_min, base10 + em->dangle5[pb][stb][st1b] + dp[piv + 1][en][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        u_min = std::min(u_min, base11 + em->terminal[pl1b][pb][stb][st1b] + right_unpaired);
        u2_min =
            std::min(u2_min, base11 + em->terminal[pl1b][pb][stb][st1b] + dp[piv + 1][en][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + em->MismatchCoaxial(pl1b, pb, stb, st1b) +
            std::min(dp[piv + 1][en][DP_U_WC], dp[piv + 1][en][DP_U_GU]);
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);

        // (   )<.(   ). > Right coax forward and backward
        val = base00 + dp[piv + 1][en][DP_U_RC];
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
        rcoax_min =
            std::min(rcoax_min, base11 + em->MismatchCoaxial(pl1b, pb, stb, st1b) + right_unpaired);

        // (   )(<   ) > Flush coax - U
        val = base01 + em->stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC];
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
        if (IsGu(pb)) {
          val = base01 + em->stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU];
          u_min = std::min(u_min, val);
          u2_min = std::min(u2_min, val);
        }
      }

      dp[st][en][DP_U] = u_min;
      dp[st][en][DP_U2] = u2_min;
      dp[st][en][DP_U_WC] = wc_min;
      dp[st][en][DP_U_GU] = gu_min;
      dp[st][en][DP_U_RC] = rcoax_min;
    }
  }
}

}  // namespace mrna::md::t22
