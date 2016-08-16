#include "fold1.h"

namespace memerna {
namespace fold {

const energy_t MULTILOOP_A = 93;

const energy_t AUGUBRANCH[4][4] = {
    {-6, -6, -6, 5 - 6},
    {-6, -6, -6, -6},
    {-6, -6, -6, 5 - 6},
    {5 - 6, -6, 5 - 6, -6}
};

energy_t FastTwoLoop(int ost, int oen, int ist, int ien) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0)
    return stacking_e[r[ost]][r[ist]][r[ien]][r[oen]];
  if (toplen == 0 || botlen == 0)
    return energy::BulgeEnergy(ost, oen, ist, ien);
  if (toplen == 1 && botlen == 1)
    return internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];

  static_assert(constants::TWOLOOP_MAX_SZ <= INITIATION_CACHE_SZ, "initiation cache not large enough");
  energy_t energy = internal_init[toplen + botlen] + std::min(std::abs(toplen - botlen) * internal_asym, constants::NINIO_MAX_ASYM);

  if (IsAuGu(r[ost], r[oen]))
    energy += internal_augu_penalty;
  if (IsAuGu(r[ist], r[ien]))
    energy += internal_augu_penalty;

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];

  return energy;
}

array3d_t<energy_t, DP_SIZE> ComputeTables1() {
  int N = int(r.size());
  // Automatically initialised to MAX_E.
  array3d_t<energy_t, DP_SIZE> arr(r.size() + 1);
  //std::vector<int> candidates;

  static_assert(constants::HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int st = N - 1; st >= 0; --st) {
    //candidates.clear();
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      base_t stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];

      // Update paired - only if can actually pair.
      if (CanPair(r[st], r[en]) && IsNotLonely(st, en)) {
        energy_t p_min = constants::MAX_E;
        int max_inter = std::min(constants::TWOLOOP_MAX_SZ, en - st - constants::HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (arr[ist][ien][DP_P] < constants::CAP_E)
              p_min = std::min(p_min, FastTwoLoop(st, en, ist, ien) + arr[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        p_min = std::min(p_min, energy::HairpinEnergy(st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = AUGUBRANCH[stb][enb] + MULTILOOP_A;

        // (<   ><   >)
        p_min = std::min(p_min, base_branch_cost + arr[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        p_min = std::min(p_min, base_branch_cost + arr[st + 2][en - 1][DP_U2] + dangle3_e[stb][st1b][enb]);
        // (<   ><   >5) 5'
        p_min = std::min(p_min, base_branch_cost + arr[st + 1][en - 2][DP_U2] + dangle5_e[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        p_min = std::min(p_min, base_branch_cost + arr[st + 2][en - 2][DP_U2] + terminal_e[stb][st1b][en1b][enb]);

        for (int piv = st + constants::HAIRPIN_MIN_SZ + 2; piv < en - constants::HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          base_t pl1b = r[piv - 1], plb = r[piv], prb = r[piv + 1], pr1b = r[piv + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          auto outer_coax = energy::MismatchMediatedCoaxialEnergy(stb, st1b, en1b, enb);
          p_min = std::min(p_min, base_branch_cost + arr[st + 2][piv][DP_P] +
              AUGUBRANCH[st2b][plb] + arr[piv + 1][en - 2][DP_U] + outer_coax);
          // (.   (   ).) Right outer coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 2][piv][DP_U] +
              AUGUBRANCH[prb][en2b] + arr[piv + 1][en - 2][DP_P] + outer_coax);

          // (.(   ).   ) Left right coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 2][piv - 1][DP_P] +
              AUGUBRANCH[st2b][pl1b] + arr[piv + 1][en - 1][DP_U] +
              energy::MismatchMediatedCoaxialEnergy(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 1][piv][DP_U] +
              AUGUBRANCH[pr1b][en2b] + arr[piv + 2][en - 2][DP_P] +
              energy::MismatchMediatedCoaxialEnergy(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 1][piv][DP_P] +
              AUGUBRANCH[st1b][plb] + arr[piv + 1][en - 1][DP_U] + stacking_e[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          p_min = std::min(p_min, base_branch_cost + arr[st + 1][piv][DP_U] +
              AUGUBRANCH[prb][en1b] + arr[piv + 1][en - 1][DP_P] + stacking_e[stb][prb][en1b][enb]);
        }

        arr[st][en][DP_P] = p_min;
      }
      energy_t u_min = constants::MAX_E, u2_min = constants::MAX_E,
          rcoax_min = constants::MAX_E, wc_min = constants::MAX_E, gu_min = constants::MAX_E;
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u_min = std::min(u_min, arr[st + 1][en][DP_U]);
        u2_min = std::min(u2_min, arr[st + 1][en][DP_U2]);
      }
      // Pair here.
      u_min = std::min(u_min, arr[st][en][DP_P] + AUGUBRANCH[stb][enb]);
      for (int piv = st + constants::HAIRPIN_MIN_SZ + 2; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        auto pb = r[piv], pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[st][piv][DP_P] + energy::AuGuPenalty(st, piv) + multiloop_hack_b;
        auto base01 = arr[st][piv - 1][DP_P] + energy::AuGuPenalty(st, piv - 1) + multiloop_hack_b;
        auto base10 = arr[st + 1][piv][DP_P] + energy::AuGuPenalty(st + 1, piv) + multiloop_hack_b;
        auto base11 = arr[st + 1][piv - 1][DP_P] + energy::AuGuPenalty(st + 1, piv - 1) + multiloop_hack_b;
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
        u_min = std::min(u_min, base01 + dangle3_e[pl1b][pb][stb] + right_unpaired);
        u2_min = std::min(u2_min, base01 + dangle3_e[pl1b][pb][stb] + arr[piv + 1][en][DP_U]);
        // 5(   )<   > 5' - U
        u_min = std::min(u_min, base10 + dangle5_e[pb][stb][st1b] + right_unpaired);
        u2_min = std::min(u2_min, base10 + dangle5_e[pb][stb][st1b] + arr[piv + 1][en][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        u_min = std::min(u_min, base11 + terminal_e[pl1b][pb][stb][st1b] + right_unpaired);
        u2_min = std::min(u2_min, base11 + terminal_e[pl1b][pb][stb][st1b] + arr[piv + 1][en][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + energy::MismatchMediatedCoaxialEnergy(pl1b, pb, stb, st1b) +
            std::min(arr[piv + 1][en][DP_U_WC], arr[piv + 1][en][DP_U_GU]);
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);

        // (   ).<(   ). > Right coax forward and backward
        val = base01 + arr[piv + 1][en][DP_U_RCOAX];
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
        if (st > 0)
          rcoax_min = std::min(rcoax_min, base01 + energy::MismatchMediatedCoaxialEnergy(
              pl1b, pb, r[st - 1], stb) + right_unpaired);

        // There has to be remaining bases to even have a chance at these cases.
        if (piv < en) {
          auto pr1b = r[piv + 1];
          // (   )<(   ) > Flush coax - U
          val = base00 + stacking_e[pb][pr1b][pr1b ^ 3][stb] + arr[piv + 1][en][DP_U_WC];
          u_min = std::min(u_min, val);
          u2_min = std::min(u2_min, val);
          if (pr1b == G || pr1b == U) {
            val = base00 + stacking_e[pb][pr1b][pr1b ^ 1][stb] + arr[piv + 1][en][DP_U_GU];
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
//      if (arr[st][en][DP_P] < arr[st][en][DP_U])
//        candidates.push_back(en);
    }
  }
  return arr;
}

}
}
