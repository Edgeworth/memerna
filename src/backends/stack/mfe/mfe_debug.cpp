// Copyright 2023 Eliot Courtney.
#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>

#include "api/energy/energy_cfg.h"
#include "backends/common/base/dp.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/mfe.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/array.h"
#include "util/error.h"

namespace mrna::md::stack {

using base::DP_P;
using base::DP_U;
using base::DP_U2;
using base::DP_U_GU;
using base::DP_U_RC;
using base::DP_U_WC;

namespace {

struct MfeInternal {
  const Primary& r;
  const Model& m;
  int N;
  base::DpArray& dp;
  Array2D<Energy>& nostack;
  Array3D<Energy>& penult;

  MfeInternal(const Primary& r_, const Model::Ptr& m_, DpState& state_)
      : r(r_), m(*m_), N(static_cast<int>(r_.size())), dp(state_.base.dp), nostack(state_.nostack),
        penult(state_.penult) {}

  void Compute() {
    static_assert(
        HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

    static thread_local const erg::EnergyCfgSupport support{
        .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
        .bulge_states{false, true},
        .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::NONE},
    };
    support.VerifySupported(funcname(), m.cfg());
    m.pf.Verify(r);

    spdlog::debug("stack {} with cfg {}", funcname(), m.cfg());

    dp = base::DpArray(r.size() + 1, MAX_E);
    nostack = Array2D(r.size() + 1, MAX_E);
    penult = Array3D(r.size() + 1, MAX_E);

    for (int st = N - 1; st >= 0; --st) {
      for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        const Base stb = r[st];
        const Base st1b = r[st + 1];
        const Base st2b = r[st + 2];
        const Base enb = r[en];
        const Base en1b = r[en - 1];
        const Base en2b = r[en - 2];

        if (m.CanPair(r, st, en)) {
          Energy stack_min = MAX_E;

          {
            const int max_stack = en - st - HAIRPIN_MIN_SZ + 1;
            // Try all stacks of each length, with or without a 1 nuc bulge loop intercedeing.
            const Energy bulge_left = m.Bulge(r, st, en, st + 2, en - 1);
            const Energy bulge_right = m.Bulge(r, st, en, st + 1, en - 2);

            // Try stems with a specific length.
            for (int length = 2; 2 * length <= max_stack; ++length) {
              // Update our DP at (st, en) - no bulge:
              if (m.CanPair(r, st + 1, en - 1)) {
                // Include the paired pseudofree energy value for all except the last pair in this
                // stack. The last pair will be handled by whatever starts after this stack.
                auto none = m.stack[r[st]][r[st + 1]][r[en - 1]][r[en]] + m.pf.Paired(st, en);
                // Try ending the stack without a bulge loop.
                if (length == 2) {
                  none += nostack[st + 1][en - 1] +
                      m.penultimate_stack[r[st]][r[st + 1]][r[en - 1]][r[en]];
                } else {
                  none += penult[st + 1][en - 1][length - 1];
                }
                penult[st][en][length] = std::min(penult[st][en][length], none);
                stack_min = std::min(stack_min, none + m.penultimate_stack[en1b][enb][stb][st1b]);
              }

              // Left bulge:
              if (m.CanPair(r, st + 2, en - 1)) {
                auto left = bulge_left;
                // Try ending the stack without a bulge loop.
                if (length == 2) {
                  left += nostack[st + 2][en - 1] +
                      m.penultimate_stack[r[st]][r[st + 2]][r[en - 1]][r[en]];
                } else {
                  left += penult[st + 2][en - 1][length - 1];
                }
                penult[st][en][length] = std::min(penult[st][en][length], left);
                stack_min = std::min(stack_min, left + m.penultimate_stack[en1b][enb][stb][st2b]);
              }

              // Right bulge:
              if (m.CanPair(r, st + 1, en - 2)) {
                auto right = bulge_right;
                if (length == 2) {
                  right += nostack[st + 1][en - 2] +
                      m.penultimate_stack[r[st]][r[st + 1]][r[en - 2]][r[en]];
                } else {
                  right += penult[st + 1][en - 2][length - 1];
                }
                penult[st][en][length] = std::min(penult[st][en][length], right);

                stack_min = std::min(stack_min, right + m.penultimate_stack[en2b][enb][stb][st1b]);
              }
            }
          }

          Energy nostack_min = MAX_E;

          const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
          for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
            for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
              // Try all internal loops. We don't check stacks or 1 nuc bulge loops.
              if (dp[ist][ien][DP_P] < CAP_E && ist - st + en - ien > 3)
                nostack_min =
                    std::min(nostack_min, m.TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P]);
            }
          }
          // Hairpin loops.
          nostack_min = std::min(nostack_min, m.Hairpin(r, st, en));

          // Multiloops. Look at range [st + 1, en - 1].
          // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
          const auto base_branch_cost = m.AuGuPenalty(stb, enb) + m.pf.Paired(st, en) +
              m.multiloop_hack_a + m.multiloop_hack_b;

          // (<   ><   >)
          nostack_min = std::min(nostack_min, base_branch_cost + dp[st + 1][en - 1][DP_U2]);

          if (m.cfg().UseDangleMismatch()) {
            // (3<   ><   >) 3'
            nostack_min = std::min(nostack_min,
                base_branch_cost + dp[st + 2][en - 1][DP_U2] + m.dangle3[stb][st1b][enb] +
                    m.pf.Unpaired(st + 1));
            // (<   ><   >5) 5'
            nostack_min = std::min(nostack_min,
                base_branch_cost + dp[st + 1][en - 2][DP_U2] + m.dangle5[stb][en1b][enb] +
                    m.pf.Unpaired(en - 1));
            // (.<   ><   >.) Terminal mismatch
            nostack_min = std::min(nostack_min,
                base_branch_cost + dp[st + 2][en - 2][DP_U2] + m.terminal[stb][st1b][en1b][enb] +
                    m.pf.Unpaired(st + 1) + m.pf.Unpaired(en - 1));
          }

          if (m.cfg().UseCoaxialStacking()) {
            const auto outer_coax = m.MismatchCoaxial(stb, st1b, en1b, enb) +
                m.pf.Unpaired(st + 1) + m.pf.Unpaired(en - 1);
            for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
              // Paired coaxial stacking cases:
              const Base pl1b = r[piv - 1];
              const Base plb = r[piv];
              const Base prb = r[piv + 1];
              const Base pr1b = r[piv + 2];
              //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
              // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

              // (.(   )   .) Left outer coax - P
              nostack_min = std::min(nostack_min,
                  base_branch_cost + dp[st + 2][piv][DP_P] + m.multiloop_hack_b +
                      m.AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax);
              // (.   (   ).) Right outer coax
              nostack_min = std::min(nostack_min,
                  base_branch_cost + dp[st + 2][piv][DP_U] + m.multiloop_hack_b +
                      m.AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax);

              // (.(   ).   ) Left inner coax
              nostack_min = std::min(nostack_min,
                  base_branch_cost + dp[st + 2][piv - 1][DP_P] + m.multiloop_hack_b +
                      m.AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
                      m.MismatchCoaxial(pl1b, plb, st1b, st2b) + m.pf.Unpaired(st + 1) +
                      m.pf.Unpaired(piv));
              // (   .(   ).) Right inner coax
              nostack_min = std::min(nostack_min,
                  base_branch_cost + dp[st + 1][piv][DP_U] + m.multiloop_hack_b +
                      m.AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
                      m.MismatchCoaxial(en2b, en1b, prb, pr1b) + m.pf.Unpaired(piv + 1) +
                      m.pf.Unpaired(en - 1));

              // ((   )   ) Left flush coax
              nostack_min = std::min(nostack_min,
                  base_branch_cost + dp[st + 1][piv][DP_P] + m.multiloop_hack_b +
                      m.AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] +
                      m.stack[stb][st1b][plb][enb]);
              // (   (   )) Right flush coax
              nostack_min = std::min(nostack_min,
                  base_branch_cost + dp[st + 1][piv][DP_U] + m.multiloop_hack_b +
                      m.AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] +
                      m.stack[stb][prb][en1b][enb]);
            }
          }

          dp[st][en][DP_P] = std::min(stack_min, nostack_min);
          nostack[st][en] = nostack_min;
        }
        Energy u_min = MAX_E;
        Energy u2_min = MAX_E;
        Energy rcoax_min = MAX_E;
        Energy wc_min = MAX_E;
        Energy gu_min = MAX_E;
        // Update unpaired.
        // Choose `st` to be unpaired.
        if (st + 1 < en) {
          u_min = std::min(u_min, dp[st + 1][en][DP_U] + m.pf.Unpaired(st));
          u2_min = std::min(u2_min, dp[st + 1][en][DP_U2] + m.pf.Unpaired(st));
        }
        for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
          //   (   .   )<   (
          // stb pl1b pb   pr1b
          const auto pb = r[piv];
          const auto pl1b = r[piv - 1];
          // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the
          // right.
          const auto base00 = dp[st][piv][DP_P] + m.multiloop_hack_b + m.AuGuPenalty(stb, pb);
          const auto base01 = dp[st][piv - 1][DP_P] + m.multiloop_hack_b + m.AuGuPenalty(stb, pl1b);
          const auto base10 = dp[st + 1][piv][DP_P] + m.multiloop_hack_b + m.AuGuPenalty(st1b, pb);
          const auto base11 =
              dp[st + 1][piv - 1][DP_P] + m.multiloop_hack_b + m.AuGuPenalty(st1b, pl1b);
          // Min is for either placing another unpaired or leaving it as nothing.
          const auto right_unpaired =
              std::min(dp[piv + 1][en][DP_U], m.pf.UnpairedCum(piv + 1, en));

          // (   )<   > - U, U_WC?, U_GU?
          u2_min = std::min(u2_min, base00 + dp[piv + 1][en][DP_U]);
          auto val = base00 + right_unpaired;
          u_min = std::min(u_min, val);
          if (IsGuPair(stb, pb))
            gu_min = std::min(gu_min, val);
          else
            wc_min = std::min(wc_min, val);

          if (m.cfg().UseDangleMismatch()) {
            // (   )3<   > 3' - U
            u_min = std::min(
                u_min, base01 + m.dangle3[pl1b][pb][stb] + m.pf.Unpaired(piv) + right_unpaired);
            u2_min = std::min(u2_min,
                base01 + m.dangle3[pl1b][pb][stb] + m.pf.Unpaired(piv) + dp[piv + 1][en][DP_U]);
            // 5(   )<   > 5' - U
            u_min = std::min(
                u_min, base10 + m.dangle5[pb][stb][st1b] + m.pf.Unpaired(st) + right_unpaired);
            u2_min = std::min(u2_min,
                base10 + m.dangle5[pb][stb][st1b] + m.pf.Unpaired(st) + dp[piv + 1][en][DP_U]);
            // .(   ).<   > Terminal mismatch - U
            u_min = std::min(u_min,
                base11 + m.terminal[pl1b][pb][stb][st1b] + m.pf.Unpaired(st) + m.pf.Unpaired(piv) +
                    right_unpaired);
            u2_min = std::min(u2_min,
                base11 + m.terminal[pl1b][pb][stb][st1b] + m.pf.Unpaired(st) + m.pf.Unpaired(piv) +
                    dp[piv + 1][en][DP_U]);
          }

          if (m.cfg().UseCoaxialStacking()) {
            // .(   ).<(   ) > Left coax - U
            val = base11 + m.MismatchCoaxial(pl1b, pb, stb, st1b) + m.pf.Unpaired(st) +
                m.pf.Unpaired(piv) + std::min(dp[piv + 1][en][DP_U_WC], dp[piv + 1][en][DP_U_GU]);
            u_min = std::min(u_min, val);
            u2_min = std::min(u2_min, val);

            // (   )<.(   ). > Right coax forward and backward
            val = base00 + dp[piv + 1][en][DP_U_RC];
            u_min = std::min(u_min, val);
            u2_min = std::min(u2_min, val);
            rcoax_min = std::min(rcoax_min,
                base11 + m.MismatchCoaxial(pl1b, pb, stb, st1b) + m.pf.Unpaired(st) +
                    m.pf.Unpaired(piv) + right_unpaired);

            // (   )(<   ) > Flush coax - U
            val = base01 + m.stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC];
            u_min = std::min(u_min, val);
            u2_min = std::min(u2_min, val);
            if (IsGu(pb)) {
              val = base01 + m.stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU];
              u_min = std::min(u_min, val);
              u2_min = std::min(u2_min, val);
            }
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
};

}  // namespace

void MfeDebug(const Primary& r, const Model::Ptr& m, DpState& state) {
  MfeInternal(r, m, state).Compute();
}

}  // namespace mrna::md::stack
