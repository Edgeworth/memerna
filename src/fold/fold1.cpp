#include "fold1.h"
#include "fold_globals.h"

namespace memerna {
namespace fold {

using constants::MAX_E;

array3d_t<energy_t, DP_SIZE> ComputeTables1() {
  InitFold();
  int N = int(r.size());
  // Automatically initialised to MAX_E.
  array3d_t<energy_t, DP_SIZE> arr(r.size() + 1);

  static_assert(constants::HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      base_t stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];

      // Update paired - only if can actually pair.
      if (CanPair(r[st], r[en]) && IsNotLonely(st, en)) {
        energy_t p_min = MAX_E;
        int max_inter = std::min(constants::TWOLOOP_MAX_SZ, en - st - constants::HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (arr[ist][ien][DP_P] < constants::CAP_E)
              p_min = std::min(p_min, FastTwoLoop(st, en, ist, ien) + arr[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        p_min = std::min(p_min, energy::Hairpin(st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = g_augubranch[stb][enb] + g_multiloop_hack_a;

        // (<   ><   >)
        p_min = std::min(p_min, base_branch_cost + arr[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        p_min = std::min(p_min, base_branch_cost + arr[st + 2][en - 1][DP_U2] + g_dangle3_e[stb][st1b][enb]);
        // (<   ><   >5) 5'
        p_min = std::min(p_min, base_branch_cost + arr[st + 1][en - 2][DP_U2] + g_dangle5_e[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        p_min = std::min(p_min, base_branch_cost + arr[st + 2][en - 2][DP_U2] + g_terminal[stb][st1b][en1b][enb]);

        for (int piv = st + constants::HAIRPIN_MIN_SZ + 2; piv < en - constants::HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          base_t pl1b = r[piv - 1], plb = r[piv], prb = r[piv + 1], pr1b = r[piv + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          auto outer_coax = energy::MismatchCoaxial(stb, st1b, en1b, enb);
          p_min = std::min(p_min, base_branch_cost + arr[st + 2][piv][DP_P] +
              g_augubranch[st2b][plb] + arr[piv + 1][en - 2][DP_U] + outer_coax);
          // (.   (   ).) Right outer coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 2][piv][DP_U] +
              g_augubranch[prb][en2b] + arr[piv + 1][en - 2][DP_P] + outer_coax);

          // (.(   ).   ) Left right coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 2][piv - 1][DP_P] +
              g_augubranch[st2b][pl1b] + arr[piv + 1][en - 1][DP_U] +
              energy::MismatchCoaxial(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 1][piv][DP_U] +
              g_augubranch[pr1b][en2b] + arr[piv + 2][en - 2][DP_P] +
              energy::MismatchCoaxial(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 1][piv][DP_P] +
              g_augubranch[st1b][plb] + arr[piv + 1][en - 1][DP_U] + g_stack[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 1][piv][DP_U] +
              g_augubranch[prb][en1b] + arr[piv + 1][en - 1][DP_P] + g_stack[stb][prb][en1b][enb]);
        }

        arr[st][en][DP_P] = p_min;
      }
      energy_t u_min = MAX_E, u2_min = MAX_E,
          rcoax_min = MAX_E, wc_min = MAX_E, gu_min = MAX_E;
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u_min = std::min(u_min, arr[st + 1][en][DP_U]);
        u2_min = std::min(u2_min, arr[st + 1][en][DP_U2]);
      }
      for (int piv = st + constants::HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        auto pb = r[piv], pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[st][piv][DP_P] + energy::AuGuPenalty(stb, pb) + g_multiloop_hack_b;
        auto base01 = arr[st][piv - 1][DP_P] + energy::AuGuPenalty(stb, pl1b) + g_multiloop_hack_b;
        auto base10 = arr[st + 1][piv][DP_P] + energy::AuGuPenalty(st1b, pb) + g_multiloop_hack_b;
        auto base11 = arr[st + 1][piv - 1][DP_P] + energy::AuGuPenalty(st1b, pl1b) + g_multiloop_hack_b;
        // Min is for either placing another unpaired or leaving it as nothing.
        auto right_unpaired = std::min(arr[piv + 1][en][DP_U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        u2_min = std::min(u2_min, base00 + arr[piv + 1][en][DP_U]);
        auto val = base00 + right_unpaired;
        u_min = std::min(u_min, val);
        if (IsGu(stb, pb))
          gu_min = std::min(gu_min, val);
        else
          wc_min = std::min(wc_min, val);

        // (   )3<   > 3' - U
        u_min = std::min(u_min, base01 + g_dangle3_e[pl1b][pb][stb] + right_unpaired);
        u2_min = std::min(u2_min, base01 + g_dangle3_e[pl1b][pb][stb] + arr[piv + 1][en][DP_U]);
        // 5(   )<   > 5' - U
        u_min = std::min(u_min, base10 + g_dangle5_e[pb][stb][st1b] + right_unpaired);
        u2_min = std::min(u2_min, base10 + g_dangle5_e[pb][stb][st1b] + arr[piv + 1][en][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        u_min = std::min(u_min, base11 + g_terminal[pl1b][pb][stb][st1b] + right_unpaired);
        u2_min = std::min(u2_min, base11 + g_terminal[pl1b][pb][stb][st1b] + arr[piv + 1][en][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + energy::MismatchCoaxial(pl1b, pb, stb, st1b) +
            std::min(arr[piv + 1][en][DP_U_WC], arr[piv + 1][en][DP_U_GU]);
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);

        // (   ).<(   ). > Right coax forward and backward
        val = base01 + arr[piv + 1][en][DP_U_RCOAX];
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
        if (st > 0)
          rcoax_min = std::min(rcoax_min, base01 + energy::MismatchCoaxial(
              pl1b, pb, r[st - 1], stb) + right_unpaired);

        // There has to be remaining bases to even have a chance at these cases.
        if (piv < en) {
          auto pr1b = r[piv + 1];
          // (   )<(   ) > Flush coax - U
          val = base00 + g_stack[pb][pr1b][pr1b ^ 3][stb] + arr[piv + 1][en][DP_U_WC];
          u_min = std::min(u_min, val);
          u2_min = std::min(u2_min, val);
          if (pr1b == G || pr1b == U) {
            val = base00 + g_stack[pb][pr1b][pr1b ^ 1][stb] + arr[piv + 1][en][DP_U_GU];
            u_min = std::min(u_min, val);
            u2_min = std::min(u2_min, val);
          }
        }
      }

      arr[st][en][DP_U] = u_min;
      arr[st][en][DP_U2] = u2_min;
      arr[st][en][DP_U_WC] = wc_min;
      arr[st][en][DP_U_GU] = gu_min;
      arr[st][en][DP_U_RCOAX] = rcoax_min;
    }
  }
  return arr;
}

}
}
