#include "fold/fold.h"
#include "array.h"

namespace memerna {
namespace fold {

#define UPDATE_CACHE(a, value) \
  do { \
    energy_t macro_upd_value_ = (value); \
    if (macro_upd_value_ < constants::CAP_E && macro_upd_value_ < arr[sz][st][a]) { \
      /*printf("Upd %d %d %d %d => %d " #value "\n", st, en, a, arr[sz][st][a], macro_upd_value_);*/ \
      arr[sz][st][a] = macro_upd_value_; \
    } \
  } while (0)


array3d_t<energy_t, DP_SIZE> ComputeTables_Slow() {
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
        for (int isz = 0; isz <= std::min(constants::TWOLOOP_MAX_SZ, sz - 4 - constants::HAIRPIN_MIN_SZ); ++isz) {
          for (int ist = st + 1; ist < st + isz + 2; ++ist) {
            int ien = en - (st + isz - ist) - 2;
            if (arr[ien - ist + 1][ist][DP_P] < constants::CAP_E)
              UPDATE_CACHE(DP_P, energy::TwoLoop(st, en, ist, ien) + arr[ien - ist + 1][ist][DP_P]);
          }
        }
        // Hairpin loops.
        UPDATE_CACHE(DP_P, energy::HairpinEnergy(st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = energy::AuGuPenalty(st, en) + multiloop_hack_a + multiloop_hack_b;

        // No stacking case.
        UPDATE_CACHE(DP_P, base_branch_cost + arr[sz - 2][st + 1][DP_U2]);
        // (3<   ><   >) 3'
        UPDATE_CACHE(DP_P, base_branch_cost + arr[sz - 3][st + 2][DP_U2] + dangle3_e[stb][st1b][enb]);
        // (<   ><   >5) 5'
        UPDATE_CACHE(DP_P, base_branch_cost + arr[sz - 3][st + 1][DP_U2] + dangle5_e[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        UPDATE_CACHE(DP_P, base_branch_cost + arr[sz - 4][st + 2][DP_U2] + terminal_e[stb][st1b][en1b][enb]);

        for (int lpivsz = constants::HAIRPIN_MIN_SZ + 2; lpivsz < sz - 3 - constants::HAIRPIN_MIN_SZ; ++lpivsz) {
          int rpivsz = sz - lpivsz - 2;
          // Paired coaxial stacking cases:
          base_t pl1b = r[st + lpivsz - 1], plb = r[st + lpivsz], prb = r[st + lpivsz + 1], pr1b = r[st + lpivsz + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          auto outer_coax = energy::MismatchMediatedCoaxialEnergy(stb, st1b, en1b, enb);
          UPDATE_CACHE(DP_P, base_branch_cost + arr[lpivsz - 1][st + 2][DP_P] + multiloop_hack_b +
              energy::AuGuPenalty(st + 2, st + lpivsz) + arr[rpivsz - 1][st + 1 + lpivsz][DP_U] + outer_coax);
          // (.   (   ).) Right outer coax
          UPDATE_CACHE(DP_P, base_branch_cost + arr[lpivsz - 1][st + 2][DP_U] + multiloop_hack_b +
              energy::AuGuPenalty(st + 1 + lpivsz, en - 2) + arr[rpivsz - 1][st + 1 + lpivsz][DP_P] + outer_coax);

          // (.(   ).   ) Left right coax
          UPDATE_CACHE(DP_P, base_branch_cost + arr[lpivsz - 2][st + 2][DP_P] + multiloop_hack_b +
              energy::AuGuPenalty(st + 2, st + lpivsz - 1) + arr[rpivsz][st + 1 + lpivsz][DP_U] +
              energy::MismatchMediatedCoaxialEnergy(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          UPDATE_CACHE(DP_P, base_branch_cost + arr[lpivsz][st + 1][DP_U] + multiloop_hack_b +
              energy::AuGuPenalty(st + 2 + lpivsz, en - 2) + arr[rpivsz - 2][st + 2 + lpivsz][DP_P] +
              energy::MismatchMediatedCoaxialEnergy(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          UPDATE_CACHE(DP_P, base_branch_cost + arr[lpivsz][st + 1][DP_P] +
              multiloop_hack_b + energy::AuGuPenalty(st + 1, st + lpivsz) +
              arr[rpivsz][st + 1 + lpivsz][DP_U] + stacking_e[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          UPDATE_CACHE(DP_P, base_branch_cost + arr[lpivsz][st + 1][DP_U] +
              multiloop_hack_b + energy::AuGuPenalty(st + 1 + lpivsz, en - 1) +
              arr[rpivsz][st + 1 + lpivsz][DP_P] + stacking_e[stb][prb][en1b][enb]);
        }
      }

      // Update unpaired.
      // Choose |st| to be unpaired.
      if (sz) {
        UPDATE_CACHE(DP_U, arr[sz - 1][st + 1][DP_U]);
        UPDATE_CACHE(DP_U2, arr[sz - 1][st + 1][DP_U2]);
      }
      // Pair here.
      UPDATE_CACHE(DP_U, arr[sz][st][DP_P] + multiloop_hack_b + energy::AuGuPenalty(st, en));
      for (int lpivsz = constants::HAIRPIN_MIN_SZ + 2; lpivsz <= sz; ++lpivsz) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        int rpivsz = sz - lpivsz;
        auto pb = r[st + lpivsz - 1], pl1b = r[st + lpivsz - 2];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[lpivsz][st][DP_P] + energy::AuGuPenalty(st, st + lpivsz - 1) + multiloop_hack_b;
        auto base01 = arr[lpivsz - 1][st][DP_P] + energy::AuGuPenalty(st, st + lpivsz - 2) + multiloop_hack_b;
        auto base10 = arr[lpivsz - 1][st + 1][DP_P] + energy::AuGuPenalty(st + 1, st + lpivsz - 1) + multiloop_hack_b;
        auto base11 = arr[lpivsz - 2][st + 1][DP_P] + energy::AuGuPenalty(st + 1, st + lpivsz - 2) + multiloop_hack_b;
        // Min is for either placing another unpaired or leaving it as nothing.
        auto right_unpaired = std::min(arr[rpivsz][st + lpivsz][DP_U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        UPDATE_CACHE(DP_U2, base00 + arr[rpivsz][st + lpivsz][DP_U]);
        auto val = base00 + right_unpaired;
        UPDATE_CACHE(DP_U, val);
        if (IsGu(stb, pb))
          UPDATE_CACHE(DP_U_GU, val);
        else
          UPDATE_CACHE(DP_U_WC, val);

        // (   )3<   > 3' - U
        UPDATE_CACHE(DP_U, base01 + dangle3_e[pl1b][pb][stb] + right_unpaired);
        UPDATE_CACHE(DP_U2, base01 + dangle3_e[pl1b][pb][stb] + arr[rpivsz][st + lpivsz][DP_U]);
        // 5(   )<   > 5' - U
        UPDATE_CACHE(DP_U, base10 + dangle5_e[pb][stb][st1b] + right_unpaired);
        UPDATE_CACHE(DP_U2, base10 + dangle5_e[pb][stb][st1b] + arr[rpivsz][st + lpivsz][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        UPDATE_CACHE(DP_U, base11 + terminal_e[pl1b][pb][stb][st1b] + right_unpaired);
        UPDATE_CACHE(DP_U2, base11 + terminal_e[pl1b][pb][stb][st1b] + arr[rpivsz][st + lpivsz][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + energy::MismatchMediatedCoaxialEnergy(pl1b, pb, stb, st1b) +
            std::min(arr[rpivsz][st + lpivsz][DP_U_WC], arr[rpivsz][st + lpivsz][DP_U_GU]);
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);

        // (   ).<(   ). > Right coax forward and backward
        val = base01 + arr[rpivsz][st + lpivsz][DP_U_RCOAX];
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);
        if (st > 0)
          UPDATE_CACHE(DP_U_RCOAX, base01 + energy::MismatchMediatedCoaxialEnergy(
              pl1b, pb, r[st - 1], stb) + right_unpaired);

        // There has to be remaining bases to even have a chance at these cases.
        if (rpivsz > 0) {
          auto pr1b = r[st + lpivsz];
          // (   )<(   ) > Flush coax - U
          val = base00 + stacking_e[pb][pr1b][pr1b ^ 3][stb] + arr[rpivsz][st + lpivsz][DP_U_WC];
          UPDATE_CACHE(DP_U, val);
          UPDATE_CACHE(DP_U2, val);
          if (pr1b == G || pr1b == U) {
            val = base00 + stacking_e[pb][pr1b][pr1b ^ 1][stb] + arr[rpivsz][st + lpivsz][DP_U_GU];
            UPDATE_CACHE(DP_U, val);
            UPDATE_CACHE(DP_U2, val);
          }
        }
      }

      en++;
    }
  }
  return arr;
}

#undef UPDATE_CACHE

}
}
