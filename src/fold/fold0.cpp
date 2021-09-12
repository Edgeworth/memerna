// Copyright 2016 E.
#include "energy/energy_globals.h"
#include "energy/fast_energy.h"
#include "fold/fold.h"

namespace mrna {
namespace fold {
namespace internal {

using energy::gem;
using energy::ViableFoldingPair;

#define UPDATE_CACHE(a, value)                                           \
  do {                                                                   \
    energy_t macro_upd_value_ = (value);                                 \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < gdp[st][en][a]) { \
      gdp[st][en][a] = macro_upd_value_;                                 \
    }                                                                    \
  } while (0)

void ComputeTables0() {
  const int N = static_cast<int>(gr.size());
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const base_t stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en],
                   en1b = gr[en - 1], en2b = gr[en - 2];

      // Update paired - only if can actually pair.
      if (ViableFoldingPair(st, en)) {
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (gdp[ist][ien][DP_P] < CAP_E)
              UPDATE_CACHE(DP_P, gem.TwoLoop(gr, st, en, ist, ien) + gdp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        UPDATE_CACHE(DP_P, gem.Hairpin(gr, st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost =
            gem.AuGuPenalty(stb, enb) + gem.multiloop_hack_a + gem.multiloop_hack_b;

        // (<   ><   >)
        UPDATE_CACHE(DP_P, base_branch_cost + gdp[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        UPDATE_CACHE(
            DP_P, base_branch_cost + gdp[st + 2][en - 1][DP_U2] + gem.dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        UPDATE_CACHE(
            DP_P, base_branch_cost + gdp[st + 1][en - 2][DP_U2] + gem.dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        UPDATE_CACHE(DP_P,
            base_branch_cost + gdp[st + 2][en - 2][DP_U2] + gem.terminal[stb][st1b][en1b][enb]);

        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          base_t pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1], pr1b = gr[piv + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);
          UPDATE_CACHE(DP_P,
              base_branch_cost + gdp[st + 2][piv][DP_P] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(st2b, plb) + gdp[piv + 1][en - 2][DP_U] + outer_coax);
          // (.   (   ).) Right outer coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + gdp[st + 2][piv][DP_U] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(prb, en2b) + gdp[piv + 1][en - 2][DP_P] + outer_coax);

          // (.(   ).   ) Left right coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + gdp[st + 2][piv - 1][DP_P] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(st2b, pl1b) + gdp[piv + 1][en - 1][DP_U] +
                  gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + gdp[st + 1][piv][DP_U] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(pr1b, en2b) + gdp[piv + 2][en - 2][DP_P] +
                  gem.MismatchCoaxial(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + gdp[st + 1][piv][DP_P] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(st1b, plb) + gdp[piv + 1][en - 1][DP_U] +
                  gem.stack[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          UPDATE_CACHE(DP_P,
              base_branch_cost + gdp[st + 1][piv][DP_U] + gem.multiloop_hack_b +
                  gem.AuGuPenalty(prb, en1b) + gdp[piv + 1][en - 1][DP_P] +
                  gem.stack[stb][prb][en1b][enb]);
        }
      }

      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        UPDATE_CACHE(DP_U, gdp[st + 1][en][DP_U]);
        UPDATE_CACHE(DP_U2, gdp[st + 1][en][DP_U2]);
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
        const auto right_unpaired = std::min(gdp[piv + 1][en][DP_U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        UPDATE_CACHE(DP_U2, base00 + gdp[piv + 1][en][DP_U]);
        auto val = base00 + right_unpaired;
        UPDATE_CACHE(DP_U, val);
        if (IsGu(stb, pb))
          UPDATE_CACHE(DP_U_GU, val);
        else
          UPDATE_CACHE(DP_U_WC, val);

        // (   )3<   > 3' - U
        UPDATE_CACHE(DP_U, base01 + gem.dangle3[pl1b][pb][stb] + right_unpaired);
        UPDATE_CACHE(DP_U2, base01 + gem.dangle3[pl1b][pb][stb] + gdp[piv + 1][en][DP_U]);
        // 5(   )<   > 5' - U
        UPDATE_CACHE(DP_U, base10 + gem.dangle5[pb][stb][st1b] + right_unpaired);
        UPDATE_CACHE(DP_U2, base10 + gem.dangle5[pb][stb][st1b] + gdp[piv + 1][en][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        UPDATE_CACHE(DP_U, base11 + gem.terminal[pl1b][pb][stb][st1b] + right_unpaired);
        UPDATE_CACHE(DP_U2, base11 + gem.terminal[pl1b][pb][stb][st1b] + gdp[piv + 1][en][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + gem.MismatchCoaxial(pl1b, pb, stb, st1b) +
            std::min(gdp[piv + 1][en][DP_U_WC], gdp[piv + 1][en][DP_U_GU]);
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);

        // (   )<.(   ). > Right coax forward and backward
        val = base00 + gdp[piv + 1][en][DP_U_RCOAX];
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);
        UPDATE_CACHE(
            DP_U_RCOAX, base11 + gem.MismatchCoaxial(pl1b, pb, stb, st1b) + right_unpaired);

        // (   )(<   ) > Flush coax - U
        val = base01 + gem.stack[pl1b][pb][pb ^ 3][stb] + gdp[piv][en][DP_U_WC];
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);
        if (pb == G || pb == U) {
          val = base01 + gem.stack[pl1b][pb][pb ^ 1][stb] + gdp[piv][en][DP_U_GU];
          UPDATE_CACHE(DP_U, val);
          UPDATE_CACHE(DP_U2, val);
        }
      }
    }
  }
}

#undef UPDATE_CACHE
}  // namespace internal
}  // namespace fold
}  // namespace mrna
