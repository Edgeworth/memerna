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
  assert(N > 0);
  // Automatically initialised to MAX_E.
  array3d_t<energy_t, DP_SIZE> arr(r.size() + 1);
  assert(arr[N - 1][0][DP_U_WC] == constants::MAX_E);

  // sz includes st and en.
  static_assert(constants::HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int sz = constants::HAIRPIN_MIN_SZ + 2; sz <= N; ++sz) {
    int en = sz - 1;
    for (int st = 0; st < N - sz + 1; ++st) {
      base_t stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];

      // Update paired - only if can actually pair.
      if (CanPair(r[st], r[en]) && IsNotLonely(st, en)) {
        energy_t p_min = constants::MAX_E;
        for (int isz = 0; isz <= std::min(constants::TWOLOOP_MAX_SZ, sz - 4 - constants::HAIRPIN_MIN_SZ); ++isz) {
          for (int ist = st + 1; ist < st + isz + 2; ++ist) {
            int ien = en - (st + isz - ist) - 2;
            if (arr[ien - ist + 1][ist][DP_P] < constants::CAP_E)
              p_min = std::min(p_min, FastTwoLoop(st, en, ist, ien) + arr[ien - ist + 1][ist][DP_P]);
          }
        }
        // Hairpin loops.
        p_min = std::min(p_min, energy::HairpinEnergy(st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = AUGUBRANCH[stb][enb] + MULTILOOP_A;

        // No stacking case.
        p_min = std::min(p_min, base_branch_cost + arr[sz - 2][st + 1][DP_U2]);
        // (3<   ><   >) 3'
        p_min = std::min(p_min, base_branch_cost + arr[sz - 3][st + 2][DP_U2] + dangle3_e[stb][st1b][enb]);
        // (<   ><   >5) 5'
        p_min = std::min(p_min, base_branch_cost + arr[sz - 3][st + 1][DP_U2] + dangle5_e[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        p_min = std::min(p_min, base_branch_cost + arr[sz - 4][st + 2][DP_U2] + terminal_e[stb][st1b][en1b][enb]);

        for (int lpivsz = constants::HAIRPIN_MIN_SZ + 2; lpivsz < sz - 3 - constants::HAIRPIN_MIN_SZ; ++lpivsz) {
          int rpivsz = sz - lpivsz - 2;
          // Paired coaxial stacking cases:
          base_t pl1b = r[st + lpivsz - 1], plb = r[st + lpivsz], prb = r[st + lpivsz + 1], pr1b = r[st + lpivsz + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          auto outer_coax = energy::MismatchMediatedCoaxialEnergy(stb, st1b, en1b, enb);
          p_min = std::min(p_min, base_branch_cost + arr[lpivsz - 1][st + 2][DP_P] +
              AUGUBRANCH[st2b][plb] + arr[rpivsz - 1][st + 1 + lpivsz][DP_U] + outer_coax);
          // (.   (   ).) Right outer coax
          p_min = std::min(p_min, base_branch_cost + arr[lpivsz - 1][st + 2][DP_U] +
              AUGUBRANCH[prb][en2b] + arr[rpivsz - 1][st + 1 + lpivsz][DP_P] + outer_coax);

          // (.(   ).   ) Left right coax
          p_min = std::min(p_min, base_branch_cost + arr[lpivsz - 2][st + 2][DP_P] +
              AUGUBRANCH[st2b][pl1b] + arr[rpivsz][st + 1 + lpivsz][DP_U] +
              energy::MismatchMediatedCoaxialEnergy(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          p_min = std::min(p_min, base_branch_cost + arr[lpivsz][st + 1][DP_U] +
              AUGUBRANCH[pr1b][en2b] + arr[rpivsz - 2][st + 2 + lpivsz][DP_P] +
              energy::MismatchMediatedCoaxialEnergy(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          p_min = std::min(p_min, base_branch_cost + arr[lpivsz][st + 1][DP_P] +
              AUGUBRANCH[st1b][plb] + arr[rpivsz][st + 1 + lpivsz][DP_U] + stacking_e[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          p_min = std::min(p_min, base_branch_cost + arr[lpivsz][st + 1][DP_U] +
              AUGUBRANCH[prb][en1b] + arr[rpivsz][st + 1 + lpivsz][DP_P] + stacking_e[stb][prb][en1b][enb]);
        }

        arr[sz][st][DP_P] = p_min;
      }
      energy_t u_min = constants::MAX_E, u2_min = constants::MAX_E,
          rcoax_min = constants::MAX_E, wc_min = constants::MAX_E, gu_min = constants::MAX_E;
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (sz) {
        u_min = std::min(u_min, arr[sz - 1][st + 1][DP_U]);
        u2_min = std::min(u2_min, arr[sz - 1][st + 1][DP_U2]);
      }
      // Pair here.
      u_min = std::min(u_min, arr[sz][st][DP_P] + AUGUBRANCH[stb][enb]);
      for (int lpivsz = constants::HAIRPIN_MIN_SZ + 2; lpivsz <= sz; ++lpivsz) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        int rpivsz = sz - lpivsz;
        auto pb = r[st + lpivsz - 1], pl1b = r[st + lpivsz - 2];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[lpivsz][st][DP_P] + AUGUBRANCH[stb][pb];
        auto base01 = arr[lpivsz - 1][st][DP_P] + AUGUBRANCH[stb][pl1b];
        auto base10 = arr[lpivsz - 1][st + 1][DP_P] + AUGUBRANCH[st1b][pb];
        auto base11 = arr[lpivsz - 2][st + 1][DP_P] + AUGUBRANCH[st1b][pl1b];
        // Min is for either placing another unpaired or leaving it as nothing.
        auto right_unpaired = std::min(arr[rpivsz][st + lpivsz][DP_U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        u2_min = std::min(u2_min, base00 + arr[rpivsz][st + lpivsz][DP_U]);
        auto val = base00 + right_unpaired;
        u_min = std::min(u_min, val);
        if (IsGu(stb, pb))
          gu_min = std::min(gu_min, val);
        else
          wc_min = std::min(wc_min, val);

        // (   )3<   > 3' - U
        u_min = std::min(u_min, base01 + dangle3_e[pl1b][pb][stb] + right_unpaired);
        u2_min = std::min(u2_min, base01 + dangle3_e[pl1b][pb][stb] + arr[rpivsz][st + lpivsz][DP_U]);
        // 5(   )<   > 5' - U
        u_min = std::min(u_min, base10 + dangle5_e[pb][stb][st1b] + right_unpaired);
        u2_min = std::min(u2_min, base10 + dangle5_e[pb][stb][st1b] + arr[rpivsz][st + lpivsz][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        u_min = std::min(u_min, base11 + terminal_e[pl1b][pb][stb][st1b] + right_unpaired);
        u2_min = std::min(u2_min, base11 + terminal_e[pl1b][pb][stb][st1b] + arr[rpivsz][st + lpivsz][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + energy::MismatchMediatedCoaxialEnergy(pl1b, pb, stb, st1b) +
            std::min(arr[rpivsz][st + lpivsz][DP_U_WC], arr[rpivsz][st + lpivsz][DP_U_GU]);
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);

        // (   ).<(   ). > Right coax forward and backward
        val = base01 + arr[rpivsz][st + lpivsz][DP_U_RCOAX];
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
        if (st > 0)
          rcoax_min = std::min(rcoax_min, base01 + energy::MismatchMediatedCoaxialEnergy(
              pl1b, pb, r[st - 1], stb) + right_unpaired);

        // There has to be remaining bases to even have a chance at these cases.
        if (rpivsz > 0) {
          auto pr1b = r[st + lpivsz];
          // (   )<(   ) > Flush coax - U
          val = base00 + stacking_e[pb][pr1b][pr1b ^ 3][stb] + arr[rpivsz][st + lpivsz][DP_U_WC];
          u_min = std::min(u_min, val);
          u2_min = std::min(u2_min, val);
          if (pr1b == G || pr1b == U) {
            val = base00 + stacking_e[pb][pr1b][pr1b ^ 1][stb] + arr[rpivsz][st + lpivsz][DP_U_GU];
            u_min = std::min(u_min, val);
            u2_min = std::min(u2_min, val);
          }
        }
      }

      arr[sz][st][DP_U] = u_min;
      arr[sz][st][DP_U2] = u2_min;
      arr[sz][st][DP_U_WC] = wc_min;
      arr[sz][st][DP_U_GU] = gu_min;
      arr[sz][st][DP_U_RCOAX] = rcoax_min;
      en++;
    }
  }
  return arr;
}

}
}
