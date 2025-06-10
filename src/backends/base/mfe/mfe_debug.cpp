// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/error.h"

namespace mrna::md::base {

#define UPDATE_CACHE(a, value)                                          \
  do {                                                                  \
    Energy macro_upd_value_ = (value);                                  \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < dp[st][en][a]) { \
      dp[st][en][a] = macro_upd_value_;                                 \
    }                                                                   \
  } while (0)

void MfeDebug(const Primary& r, const Model::Ptr& m, DpState& state) {
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false, true},
      .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::D2,
          erg::EnergyCfg::Ctd::NONE},
  };
  support.VerifySupported(funcname(), m->cfg());
  m->pf.Verify(r);

  spdlog::debug("base {} with cfg {}", funcname(), m->cfg());

  const int N = static_cast<int>(r.size());
  state.dp = DpArray(r.size() + 1, MAX_E);
  auto& dp = state.dp;

  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const Base stb = r[st];
      const Base st1b = r[st + 1];
      const Base st2b = r[st + 2];
      const Base enb = r[en];
      const Base en1b = r[en - 1];
      const Base en2b = r[en - 2];

      if (m->CanPair(r, st, en)) {
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (dp[ist][ien][DP_P] < CAP_E)
              UPDATE_CACHE(DP_P, m->TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        UPDATE_CACHE(DP_P, m->Hairpin(r, st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost =
            m->AuGuPenalty(stb, enb) + m->pf.Paired(st, en) + m->multiloop_a + m->multiloop_b;

        // (<   ><   >)
        auto val = base_branch_cost + dp[st + 1][en - 1][DP_U2];
        if (m->cfg().UseD2()) {
          // D2 can overlap terminal mismatches with anything.
          // (<   ><   >) Terminal mismatch
          val += m->terminal[stb][st1b][en1b][enb];
        }
        UPDATE_CACHE(DP_P, val);

        if (m->cfg().UseDangleMismatch()) {
          // (3<   ><   >) 3'
          UPDATE_CACHE(DP_P,
              base_branch_cost + dp[st + 2][en - 1][DP_U2] + m->dangle3[stb][st1b][enb] +
                  m->pf.Unpaired(st + 1) + m->multiloop_c);
          // (<   ><   >5) 5'
          UPDATE_CACHE(DP_P,
              base_branch_cost + dp[st + 1][en - 2][DP_U2] + m->dangle5[stb][en1b][enb] +
                  m->pf.Unpaired(en - 1) + m->multiloop_c);
          // (.<   ><   >.) Terminal mismatch
          UPDATE_CACHE(DP_P,
              base_branch_cost + dp[st + 2][en - 2][DP_U2] + m->terminal[stb][st1b][en1b][enb] +
                  m->pf.Unpaired(st + 1) + m->pf.Unpaired(en - 1) + 2 * m->multiloop_c);
        }

        if (m->cfg().UseCoaxialStacking()) {
          const auto outer_coax = m->MismatchCoaxial(stb, st1b, en1b, enb) +
              m->pf.Unpaired(st + 1) + m->pf.Unpaired(en - 1) + 2 * m->multiloop_c;
          for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
            // Paired coaxial stacking cases:
            const Base pl1b = r[piv - 1];
            const Base plb = r[piv];
            const Base prb = r[piv + 1];
            const Base pr1b = r[piv + 2];
            //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
            // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

            // (.(   )   .) Left outer coax - P
            UPDATE_CACHE(DP_P,
                base_branch_cost + dp[st + 2][piv][DP_P] + m->multiloop_b +
                    m->AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax);
            // (.   (   ).) Right outer coax
            UPDATE_CACHE(DP_P,
                base_branch_cost + dp[st + 2][piv][DP_U] + m->multiloop_b +
                    m->AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax);

            // (.(   ).   ) Left inner coax
            UPDATE_CACHE(DP_P,
                base_branch_cost + dp[st + 2][piv - 1][DP_P] + m->multiloop_b +
                    m->AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
                    m->MismatchCoaxial(pl1b, plb, st1b, st2b) + m->pf.Unpaired(st + 1) +
                    m->pf.Unpaired(piv) + 2 * m->multiloop_c);
            // (   .(   ).) Right inner coax
            UPDATE_CACHE(DP_P,
                base_branch_cost + dp[st + 1][piv][DP_U] + m->multiloop_b +
                    m->AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
                    m->MismatchCoaxial(en2b, en1b, prb, pr1b) + m->pf.Unpaired(piv + 1) +
                    m->pf.Unpaired(en - 1) + 2 * m->multiloop_c);

            // ((   )   ) Left flush coax
            UPDATE_CACHE(DP_P,
                base_branch_cost + dp[st + 1][piv][DP_P] + m->multiloop_b +
                    m->AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] +
                    m->stack[stb][st1b][plb][enb]);
            // (   (   )) Right flush coax
            UPDATE_CACHE(DP_P,
                base_branch_cost + dp[st + 1][piv][DP_U] + m->multiloop_b +
                    m->AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] +
                    m->stack[stb][prb][en1b][enb]);
          }
        }
      }

      // Update unpaired.
      // Choose `st` to be unpaired.
      if (st + 1 < en) {
        UPDATE_CACHE(DP_U, dp[st + 1][en][DP_U] + m->pf.Unpaired(st));
        UPDATE_CACHE(DP_U2, dp[st + 1][en][DP_U2] + m->pf.Unpaired(st));
      }
      // Pair here.
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = r[piv];
        const auto pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const auto base00 = dp[st][piv][DP_P] + m->AuGuPenalty(stb, pb) + m->multiloop_b;
        const auto base01 = dp[st][piv - 1][DP_P] + m->AuGuPenalty(stb, pl1b) + m->multiloop_b;
        const auto base10 = dp[st + 1][piv][DP_P] + m->AuGuPenalty(st1b, pb) + m->multiloop_b;
        const auto base11 = dp[st + 1][piv - 1][DP_P] + m->AuGuPenalty(st1b, pl1b) + m->multiloop_b;
        // Min is for either placing another unpaired or leaving it as nothing.
        const auto right_unpaired = std::min(
            dp[piv + 1][en][DP_U], m->pf.UnpairedCum(piv + 1, en) + (en - piv) * m->multiloop_c);

        // (   )<   > - U, U_WC?, U_GU?
        auto u2_val = base00 + dp[piv + 1][en][DP_U];
        auto val = base00 + right_unpaired;

        if (m->cfg().UseD2()) {
          // Note that D2 can overlap with anything.
          if (st != 0 && piv != N - 1) {
            // (   )<   > Terminal mismatch - U
            val += m->terminal[pb][r[piv + 1]][r[st - 1]][stb];
            u2_val += m->terminal[pb][r[piv + 1]][r[st - 1]][stb];
          } else if (piv != N - 1) {
            // (   )<3   > 3' - U
            val += m->dangle3[pb][r[piv + 1]][stb];
            u2_val += m->dangle3[pb][r[piv + 1]][stb];
          } else if (st != 0) {
            // 5(   )<   > 5' - U
            val += m->dangle5[pb][r[st - 1]][stb];
            u2_val += m->dangle5[pb][r[st - 1]][stb];
          }
        }

        UPDATE_CACHE(DP_U2, u2_val);
        UPDATE_CACHE(DP_U, val);
        if (IsGuPair(stb, pb))
          UPDATE_CACHE(DP_U_GU, val);
        else
          UPDATE_CACHE(DP_U_WC, val);

        if (m->cfg().UseDangleMismatch()) {
          // (   )3<   > 3' - U
          UPDATE_CACHE(DP_U,
              base01 + m->dangle3[pl1b][pb][stb] + m->pf.Unpaired(piv) + m->multiloop_c +
                  right_unpaired);
          UPDATE_CACHE(DP_U2,
              base01 + m->dangle3[pl1b][pb][stb] + m->pf.Unpaired(piv) + m->multiloop_c +
                  dp[piv + 1][en][DP_U]);
          // 5(   )<   > 5' - U
          UPDATE_CACHE(DP_U,
              base10 + m->dangle5[pb][stb][st1b] + m->pf.Unpaired(st) + m->multiloop_c +
                  right_unpaired);
          UPDATE_CACHE(DP_U2,
              base10 + m->dangle5[pb][stb][st1b] + m->pf.Unpaired(st) + m->multiloop_c +
                  dp[piv + 1][en][DP_U]);
          // .(   ).<   > Terminal mismatch - U
          UPDATE_CACHE(DP_U,
              base11 + m->terminal[pl1b][pb][stb][st1b] + m->pf.Unpaired(st) + m->pf.Unpaired(piv) +
                  2 * m->multiloop_c + right_unpaired);
          UPDATE_CACHE(DP_U2,
              base11 + m->terminal[pl1b][pb][stb][st1b] + m->pf.Unpaired(st) + m->pf.Unpaired(piv) +
                  2 * m->multiloop_c + dp[piv + 1][en][DP_U]);
        }

        if (m->cfg().UseCoaxialStacking()) {
          // .(   ).<(   ) > Left coax - U
          val = base11 + m->MismatchCoaxial(pl1b, pb, stb, st1b) + m->pf.Unpaired(st) +
              m->pf.Unpaired(piv) + 2 * m->multiloop_c +
              std::min(dp[piv + 1][en][DP_U_WC], dp[piv + 1][en][DP_U_GU]);
          UPDATE_CACHE(DP_U, val);
          UPDATE_CACHE(DP_U2, val);

          // (   )<.(   ). > Right coax forward and backward
          val = base00 + dp[piv + 1][en][DP_U_RC];
          UPDATE_CACHE(DP_U, val);
          UPDATE_CACHE(DP_U2, val);
          UPDATE_CACHE(DP_U_RC,
              base11 + m->MismatchCoaxial(pl1b, pb, stb, st1b) + m->pf.Unpaired(st) +
                  m->pf.Unpaired(piv) + 2 * m->multiloop_c + right_unpaired);

          // (   )(<   ) > Flush coax - U
          val = base01 + m->stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC];
          UPDATE_CACHE(DP_U, val);
          UPDATE_CACHE(DP_U2, val);
          if (IsGu(pb)) {
            val = base01 + m->stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU];
            UPDATE_CACHE(DP_U, val);
            UPDATE_CACHE(DP_U2, val);
          }
        }
      }
    }
  }
}

#undef UPDATE_CACHE

}  // namespace mrna::md::base
