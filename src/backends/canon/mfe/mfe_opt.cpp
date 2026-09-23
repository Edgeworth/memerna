// Copyright 2026 Eliot Courtney.
#include "backends/canon/mfe/mfe_opt.h"

#include <fmt/format.h>

#include <algorithm>
#include <string>

#include "api/energy/energy_cfg.h"
#include "backends/canon/dp.h"
#include "backends/canon/energy/model.h"
#include "backends/canon/energy/precomp.h"
#include "backends/common/static.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/error.h"
#include "util/log.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
bool MfeOpt<CFG, PseudofreeCfgT>::IsSupported(const PseudofreeCfgT& /*pf*/, std::string* reason) {
  if constexpr (CFG.lonely_pairs != erg::EnergyCfg::LonelyPairs::HEURISTIC &&
      CFG.lonely_pairs != erg::EnergyCfg::LonelyPairs::ON) {
    if (reason) *reason = fmt::format("lonely_pairs={} not supported", CFG.lonely_pairs);
    return false;
  }
  return true;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
void MfeOpt<CFG, PseudofreeCfgT>::Run(
    const Primary& r, const Model::Ptr& m, DpState& state, const PseudofreeCfgT& pf) {
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  std::string reason;
  verify(IsSupported(pf, &reason), "{} does not support the given configuration: {}", funcname(),
      reason);
  pf.Verify(r);

  logdebug("canon {} with {}, {}", funcname(), CFG, pf);

  const int N = r.size();
  const Precomp<CFG, PseudofreeCfgT> pc(Primary(r), m, pf);
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

      if (Model::CanPair<CFG>(r, st, en)) {
        Energy p_min = MAX_E;
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (dp[ist][ien][DP_P] < CAP_E)
              p_min = std::min(p_min, pc.TwoLoop(st, en, ist, ien) + dp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        p_min = std::min(p_min, pc.Hairpin(st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop
        // helix.
        const auto base_branch_cost = pc.augubranch[stb][enb] + pf.Paired(st, en) + m->multiloop_a;

        // (<   ><   >)
        auto val = base_branch_cost + dp[st + 1][en - 1][DP_U2];
        if constexpr (CFG.UseD2()) {
          // D2 can overlap terminal mismatches with anything.
          // (<   ><   >) Terminal mismatch
          val += m->terminal[stb][st1b][en1b][enb];
        }
        p_min = std::min(p_min, val);

        if constexpr (CFG.UseDangleMismatch()) {
          // (3<   ><   >) 3'
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 2][en - 1][DP_U2] + m->dangle3[stb][st1b][enb] +
                  pf.Unpaired(st + 1) + m->multiloop_c);
          // (<   ><   >5) 5'
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 1][en - 2][DP_U2] + m->dangle5[stb][en1b][enb] +
                  pf.Unpaired(en - 1) + m->multiloop_c);
          // (.<   ><   >.) Terminal mismatch
          p_min = std::min(p_min,
              base_branch_cost + dp[st + 2][en - 2][DP_U2] + m->terminal[stb][st1b][en1b][enb] +
                  pf.Unpaired(st + 1) + pf.Unpaired(en - 1) + 2 * m->multiloop_c);
        }

        if constexpr (CFG.UseCoaxialStacking()) {
          const auto outer_coax = m->MismatchCoaxial(stb, st1b, en1b, enb) + pf.Unpaired(st + 1) +
              pf.Unpaired(en - 1) + 2 * m->multiloop_c;
          for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
            // Paired coaxial stacking cases:
            const Base pl1b = r[piv - 1];
            const Base plb = r[piv];
            const Base prb = r[piv + 1];
            const Base pr1b = r[piv + 2];

            // (.(   )   .) Left outer coax - P
            p_min = std::min(p_min,
                base_branch_cost + dp[st + 2][piv][DP_P] + pc.augubranch[st2b][plb] +
                    dp[piv + 1][en - 2][DP_U] + outer_coax);
            // (.   (   ).) Right outer coax
            p_min = std::min(p_min,
                base_branch_cost + dp[st + 2][piv][DP_U] + pc.augubranch[prb][en2b] +
                    dp[piv + 1][en - 2][DP_P] + outer_coax);

            // (.(   ).   ) Left inner coax
            p_min = std::min(p_min,
                base_branch_cost + dp[st + 2][piv - 1][DP_P] + pc.augubranch[st2b][pl1b] +
                    dp[piv + 1][en - 1][DP_U] + m->MismatchCoaxial(pl1b, plb, st1b, st2b) +
                    pf.Unpaired(st + 1) + pf.Unpaired(piv) + 2 * m->multiloop_c);
            // (   .(   ).) Right inner coax
            p_min = std::min(p_min,
                base_branch_cost + dp[st + 1][piv][DP_U] + pc.augubranch[pr1b][en2b] +
                    dp[piv + 2][en - 2][DP_P] + m->MismatchCoaxial(en2b, en1b, prb, pr1b) +
                    pf.Unpaired(piv + 1) + pf.Unpaired(en - 1) + 2 * m->multiloop_c);

            // ((   )   ) Left flush coax
            p_min = std::min(p_min,
                base_branch_cost + dp[st + 1][piv][DP_P] + pc.augubranch[st1b][plb] +
                    dp[piv + 1][en - 1][DP_U] + m->stack[stb][st1b][plb][enb]);
            // (   (   )) Right flush coax
            p_min = std::min(p_min,
                base_branch_cost + dp[st + 1][piv][DP_U] + pc.augubranch[prb][en1b] +
                    dp[piv + 1][en - 1][DP_P] + m->stack[stb][prb][en1b][enb]);
          }
        }

        dp[st][en][DP_P] = p_min;
      }

      Energy cx_min = MAX_E;
      if constexpr (CFG.UseCoaxialStacking()) {
        // (   )(<   ) Flush coax - CX
        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv <= en - HAIRPIN_MIN_SZ - 1; ++piv) {
          const Base pl1b = r[piv - 1];
          const Base pb = r[piv];
          cx_min = std::min(cx_min,
              dp[st][piv - 1][DP_P] + pc.augubranch[stb][pl1b] + dp[piv][en][DP_P] +
                  pc.augubranch[pb][enb] + m->stack[pl1b][pb][enb][stb]);
        }
        // .(   ).<(   ) Left coax - CX
        for (int piv = st + HAIRPIN_MIN_SZ + 3; piv <= en - HAIRPIN_MIN_SZ - 2; ++piv) {
          const Base pl1b = r[piv - 1];
          const Base pb = r[piv];
          const Base prb = r[piv + 1];
          cx_min = std::min(cx_min,
              dp[st + 1][piv - 1][DP_P] + pc.augubranch[st1b][pl1b] + dp[piv + 1][en][DP_P] +
                  pc.augubranch[prb][enb] + m->MismatchCoaxial(pl1b, pb, stb, st1b) +
                  pf.Unpaired(st) + pf.Unpaired(piv) + 2 * m->multiloop_c);
        }
        // (   )<.(   ). Right coax - CX
        for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en - HAIRPIN_MIN_SZ - 4; ++piv) {
          const Base pb = r[piv];
          const Base prb = r[piv + 1];
          const Base pr1b = r[piv + 2];
          cx_min = std::min(cx_min,
              dp[st][piv][DP_P] + pc.augubranch[stb][pb] + dp[piv + 2][en - 1][DP_P] +
                  pc.augubranch[pr1b][en1b] + m->MismatchCoaxial(en1b, enb, prb, pr1b) +
                  pf.Unpaired(piv + 1) + pf.Unpaired(en) + 2 * m->multiloop_c);
        }
      }
      dp[st][en][DP_CX] = cx_min;

      Energy u_min = MAX_E;
      Energy u2_min = MAX_E;
      if (st + 1 < en) {
        u_min = std::min(u_min, dp[st + 1][en][DP_U] + pf.Unpaired(st) + m->multiloop_c);
        u2_min = std::min(u2_min, dp[st + 1][en][DP_U2] + pf.Unpaired(st) + m->multiloop_c);
      }

      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = r[piv];
        const auto pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on
        // the right.
        const auto base00 = dp[st][piv][DP_P] + pc.augubranch[stb][pb];
        const auto base01 = dp[st][piv - 1][DP_P] + pc.augubranch[stb][pl1b];
        const auto base10 = dp[st + 1][piv][DP_P] + pc.augubranch[st1b][pb];
        const auto base11 = dp[st + 1][piv - 1][DP_P] + pc.augubranch[st1b][pl1b];
        const auto right_unpaired = std::min(
            dp[piv + 1][en][DP_U], pf.UnpairedSum(piv + 1, en) + (en - piv) * m->multiloop_c);

        // (   )<   > - U, U2
        auto u2_val = base00 + dp[piv + 1][en][DP_U];
        auto val = base00 + right_unpaired;

        if constexpr (CFG.UseD2()) {
          // Note that D2 can overlap with anything.
          if (st != 0 && piv != N - 1) {
            // (   )<   > Terminal mismatch - U, U2
            val += m->terminal[pb][r[piv + 1]][r[st - 1]][stb];
            u2_val += m->terminal[pb][r[piv + 1]][r[st - 1]][stb];
          } else if (piv != N - 1) {
            // (   )<3   > 3' - U, U2
            val += m->dangle3[pb][r[piv + 1]][stb];
            u2_val += m->dangle3[pb][r[piv + 1]][stb];
          } else if (st != 0) {
            // 5(   )<   > 5' - U, U2
            val += m->dangle5[pb][r[st - 1]][stb];
            u2_val += m->dangle5[pb][r[st - 1]][stb];
          }
        }

        u2_min = std::min(u2_min, u2_val);
        u_min = std::min(u_min, val);

        if constexpr (CFG.UseDangleMismatch()) {
          // (   )3<   > 3' - U, U2
          u_min = std::min(u_min,
              base01 + m->dangle3[pl1b][pb][stb] + pf.Unpaired(piv) + m->multiloop_c +
                  right_unpaired);
          u2_min = std::min(u2_min,
              base01 + m->dangle3[pl1b][pb][stb] + pf.Unpaired(piv) + m->multiloop_c +
                  dp[piv + 1][en][DP_U]);
          // 5(   )<   > 5' - U, U2
          u_min = std::min(u_min,
              base10 + m->dangle5[pb][stb][st1b] + pf.Unpaired(st) + m->multiloop_c +
                  right_unpaired);
          u2_min = std::min(u2_min,
              base10 + m->dangle5[pb][stb][st1b] + pf.Unpaired(st) + m->multiloop_c +
                  dp[piv + 1][en][DP_U]);
          // .(   ).<   > Terminal mismatch - U, U2
          u_min = std::min(u_min,
              base11 + m->terminal[pl1b][pb][stb][st1b] + pf.Unpaired(st) + pf.Unpaired(piv) +
                  2 * m->multiloop_c + right_unpaired);
          u2_min = std::min(u2_min,
              base11 + m->terminal[pl1b][pb][stb][st1b] + pf.Unpaired(st) + pf.Unpaired(piv) +
                  2 * m->multiloop_c + dp[piv + 1][en][DP_U]);
        }
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        for (int piv = st + 2 * HAIRPIN_MIN_SZ + 3; piv <= en; ++piv) {
          // (   )(<   )<   > - CX
          const auto right_unpaired = std::min(
              dp[piv + 1][en][DP_U], pf.UnpairedSum(piv + 1, en) + (en - piv) * m->multiloop_c);
          const auto val = dp[st][piv][DP_CX] + right_unpaired;
          u_min = std::min(u_min, val);
          u2_min = std::min(u2_min, val);
        }
      }

      dp[st][en][DP_U] = u_min;
      dp[st][en][DP_U2] = u2_min;
    }
  }
}

#define MRNA_INSTANTIATE(CFG, PF) template class MfeOpt<CFG, PF>;
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::canon
