// Copyright 2016 Eliot Courtney.
#include "backends/baseopt/trace/trace.h"

#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "backends/common/base/dp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/error.h"

namespace mrna::md::base::opt {

TraceResult Traceback(
    const Primary& r, const Model::Ptr& m, const trace::TraceCfg& cfg, const DpState& state) {
  const int N = static_cast<int>(r.size());

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false, true},
      .ctd{erg::EnergyCfg::Ctd::ALL},
  };
  support.VerifySupported(funcname(), m->cfg());
  verify(!cfg.random, "random traceback is not supported in this energy model");

  spdlog::debug("baseopt {} with cfg {}", funcname(), m->cfg());

  const auto& [dp, ext] = state;
  TraceResult res((Secondary(N)), Ctds(N));
  std::vector<DpIndex> q;
  q.emplace_back(0, -1, EXT);
  while (!q.empty()) {
    const int st = q.back().st;
    int en = q.back().en;
    const int a = q.back().a;
    q.pop_back();

    if (en == -1) {
      // Case: No pair starting here
      if (a == EXT && st + 1 < N && ext[st + 1][EXT] == ext[st][EXT]) {
        q.emplace_back(st + 1, -1, EXT);
        goto loopend;
      }
      for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        // .   .   .   (   .   .   .   )   <   >
        //           stb  st1b   en1b  enb   rem
        const auto stb = r[st];
        const auto st1b = r[st + 1];
        const auto enb = r[en];
        const auto en1b = r[en - 1];
        const auto base00 = dp[st][en][DP_P] + m->AuGuPenalty(stb, enb);
        const auto base01 = dp[st][en - 1][DP_P] + m->AuGuPenalty(stb, en1b);
        const auto base10 = dp[st + 1][en][DP_P] + m->AuGuPenalty(st1b, enb);
        const auto base11 = dp[st + 1][en - 1][DP_P] + m->AuGuPenalty(st1b, en1b);

        // (   )<.( * ). > Right coax backward
        if (a == EXT_RC) {
          // Don't set CTDs here since they will have already been set.
          if (base11 + m->MismatchCoaxial(en1b, enb, stb, st1b) + ext[en + 1][EXT] ==
              ext[st][EXT_RC]) {
            q.emplace_back(st + 1, en - 1, DP_P);
            q.emplace_back(en + 1, -1, EXT);
            goto loopend;
          }
        }

        // EXT_RC is only for the above case.
        if (a == EXT_RC) continue;

        // (   )<   >
        auto val = base00 + ext[en + 1][EXT];
        if (val == ext[st][a] && (a != EXT_WC || IsWcPair(stb, enb)) &&
            (a != EXT_GU || IsGuPair(stb, enb))) {
          // EXT_WC and EXT_GU will have already had their ctds set.
          if (a == EXT) res.ctd[st] = CTD_UNUSED;
          q.emplace_back(st, en, DP_P);
          q.emplace_back(en + 1, -1, EXT);
          goto loopend;
        }

        // Only look at EXT from here on.
        if (a != EXT) continue;

        // (   )3<   > 3'
        if (base01 + m->dangle3[en1b][enb][stb] + ext[en + 1][EXT] == ext[st][EXT]) {
          res.ctd[st] = CTD_3_DANGLE;
          q.emplace_back(st, en - 1, DP_P);
          q.emplace_back(en + 1, -1, EXT);
          goto loopend;
        }
        // 5(   )<   > 5'
        if (base10 + m->dangle5[enb][stb][st1b] + ext[en + 1][EXT] == ext[st][EXT]) {
          res.ctd[st + 1] = CTD_5_DANGLE;
          q.emplace_back(st + 1, en, DP_P);
          q.emplace_back(en + 1, -1, EXT);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch
        if (base11 + m->terminal[en1b][enb][stb][st1b] + ext[en + 1][EXT] == ext[st][EXT]) {
          res.ctd[st + 1] = CTD_MISMATCH;
          q.emplace_back(st + 1, en - 1, DP_P);
          q.emplace_back(en + 1, -1, EXT);
          goto loopend;
        }

        if (en < N - 1) {
          // .(   ).<(   ) > Left coax  x
          val = base11 + m->MismatchCoaxial(en1b, enb, stb, st1b);
          if (val + ext[en + 1][EXT_WC] == ext[st][EXT]) {
            res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
            res.ctd[en + 1] = CTD_LCOAX_WITH_PREV;
            q.emplace_back(st + 1, en - 1, DP_P);
            q.emplace_back(en + 1, -1, EXT_WC);
            goto loopend;
          }
          if (val + ext[en + 1][EXT_GU] == ext[st][EXT]) {
            res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
            res.ctd[en + 1] = CTD_LCOAX_WITH_PREV;
            q.emplace_back(st + 1, en - 1, DP_P);
            q.emplace_back(en + 1, -1, EXT_GU);
            goto loopend;
          }

          // (   )<.(   ). > Right coax forward
          if (base00 + ext[en + 1][EXT_RC] == ext[st][EXT]) {
            res.ctd[st] = CTD_RC_WITH_NEXT;
            res.ctd[en + 2] = CTD_RC_WITH_PREV;
            q.emplace_back(st, en, DP_P);
            q.emplace_back(en + 1, -1, EXT_RC);
            goto loopend;
          }

          // (   )(<   ) > Flush coax
          if (base01 + m->stack[en1b][enb][WcPair(enb)][stb] + ext[en][EXT_WC] == ext[st][EXT]) {
            res.ctd[st] = CTD_FCOAX_WITH_NEXT;
            res.ctd[en] = CTD_FCOAX_WITH_PREV;
            q.emplace_back(st, en - 1, DP_P);
            q.emplace_back(en, -1, EXT_WC);
            goto loopend;
          }
          if (IsGu(enb) &&
              base01 + m->stack[en1b][enb][GuPair(enb)][stb] + ext[en][EXT_GU] == ext[st][EXT]) {
            res.ctd[st] = CTD_FCOAX_WITH_NEXT;
            res.ctd[en] = CTD_FCOAX_WITH_PREV;
            q.emplace_back(st, en - 1, DP_P);
            q.emplace_back(en, -1, EXT_GU);
            goto loopend;
          }
        }
      }
    } else {
      const auto stb = r[st];
      const auto st1b = r[st + 1];
      const auto st2b = r[st + 2];
      const auto enb = r[en];
      const auto en1b = r[en - 1];
      const auto en2b = r[en - 2];
      if (a == DP_P) {
        // It's paired, so add it to the folding.
        res.s[st] = en;
        res.s[en] = st;

        // Following largely matches the above DP so look up there for comments.
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (dp[ist][ien][DP_P] < CAP_E) {
              const auto val = m->TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P];
              if (val == dp[st][en][DP_P]) {
                q.emplace_back(ist, ien, DP_P);
                goto loopend;
              }
            }
          }
        }

        if (m->Hairpin(r, st, en) == dp[st][en][DP_P]) {
          goto loopend;
        }

        const auto base_branch_cost = m->AuGuPenalty(stb, enb) + m->multiloop_a + m->multiloop_b;
        // (<   ><    >)
        if (base_branch_cost + dp[st + 1][en - 1][DP_U2] == dp[st][en][DP_P]) {
          res.ctd[en] = CTD_UNUSED;
          q.emplace_back(st + 1, en - 1, DP_U2);
          goto loopend;
        }
        // (3<   ><   >) 3'
        if (base_branch_cost + dp[st + 2][en - 1][DP_U2] + m->dangle3[stb][st1b][enb] ==
            dp[st][en][DP_P]) {
          res.ctd[en] = CTD_3_DANGLE;
          q.emplace_back(st + 2, en - 1, DP_U2);
          goto loopend;
        }
        // (<   ><   >5) 5'
        if (base_branch_cost + dp[st + 1][en - 2][DP_U2] + m->dangle5[stb][en1b][enb] ==
            dp[st][en][DP_P]) {
          res.ctd[en] = CTD_5_DANGLE;
          q.emplace_back(st + 1, en - 2, DP_U2);
          goto loopend;
        }
        // (.<   ><   >.) Terminal mismatch
        if (base_branch_cost + dp[st + 2][en - 2][DP_U2] + m->terminal[stb][st1b][en1b][enb] ==
            dp[st][en][DP_P]) {
          res.ctd[en] = CTD_MISMATCH;
          q.emplace_back(st + 2, en - 2, DP_U2);
          goto loopend;
        }

        const auto outer_coax = m->MismatchCoaxial(stb, st1b, en1b, enb);
        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          const Base pl1b = r[piv - 1];
          const Base plb = r[piv];
          const Base prb = r[piv + 1];
          const Base pr1b = r[piv + 2];

          // (.(   )   .) Left outer coax - P
          if (base_branch_cost + dp[st + 2][piv][DP_P] + m->multiloop_b +
                  m->AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax ==
              dp[st][en][DP_P]) {
            res.ctd[en] = CTD_LCOAX_WITH_NEXT;
            res.ctd[st + 2] = CTD_LCOAX_WITH_PREV;
            q.emplace_back(st + 2, piv, DP_P);
            q.emplace_back(piv + 1, en - 2, DP_U);
            goto loopend;
          }
          // (.   (   ).) Right outer coax
          if (base_branch_cost + dp[st + 2][piv][DP_U] + m->multiloop_b +
                  m->AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax ==
              dp[st][en][DP_P]) {
            res.ctd[en] = CTD_RC_WITH_PREV;
            res.ctd[piv + 1] = CTD_RC_WITH_NEXT;
            q.emplace_back(st + 2, piv, DP_U);
            q.emplace_back(piv + 1, en - 2, DP_P);
            goto loopend;
          }

          // (.(   ).   ) Left inner coax
          if (base_branch_cost + dp[st + 2][piv - 1][DP_P] + m->multiloop_b +
                  m->AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
                  m->MismatchCoaxial(pl1b, plb, st1b, st2b) ==
              dp[st][en][DP_P]) {
            res.ctd[en] = CTD_RC_WITH_NEXT;
            res.ctd[st + 2] = CTD_RC_WITH_PREV;
            q.emplace_back(st + 2, piv - 1, DP_P);
            q.emplace_back(piv + 1, en - 1, DP_U);
            goto loopend;
          }
          // (   .(   ).) Right inner coax
          if (base_branch_cost + dp[st + 1][piv][DP_U] + m->multiloop_b +
                  m->AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
                  m->MismatchCoaxial(en2b, en1b, prb, pr1b) ==
              dp[st][en][DP_P]) {
            res.ctd[en] = CTD_LCOAX_WITH_PREV;
            res.ctd[piv + 2] = CTD_LCOAX_WITH_NEXT;
            q.emplace_back(st + 1, piv, DP_U);
            q.emplace_back(piv + 2, en - 2, DP_P);
            goto loopend;
          }

          // ((   )   ) Left flush coax
          if (base_branch_cost + dp[st + 1][piv][DP_P] + m->multiloop_b +
                  m->AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] +
                  m->stack[stb][st1b][plb][enb] ==
              dp[st][en][DP_P]) {
            res.ctd[en] = CTD_FCOAX_WITH_NEXT;
            res.ctd[st + 1] = CTD_FCOAX_WITH_PREV;
            q.emplace_back(st + 1, piv, DP_P);
            q.emplace_back(piv + 1, en - 1, DP_U);
            goto loopend;
          }
          // (   (   )) Right flush coax
          if (base_branch_cost + dp[st + 1][piv][DP_U] + m->multiloop_b +
                  m->AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] +
                  m->stack[stb][prb][en1b][enb] ==
              dp[st][en][DP_P]) {
            res.ctd[en] = CTD_FCOAX_WITH_PREV;
            res.ctd[piv + 1] = CTD_FCOAX_WITH_NEXT;
            q.emplace_back(st + 1, piv, DP_U);
            q.emplace_back(piv + 1, en - 1, DP_P);
            goto loopend;
          }
        }

        // Done with paired. We might not have jumped to loopend if this was a hairpin.
        continue;
      }

      // Deal with the rest of the cases:
      // Left unpaired. Either DP_U or DP_U2.
      if (st + 1 < en && (a == DP_U || a == DP_U2) && dp[st + 1][en][a] == dp[st][en][a]) {
        q.emplace_back(st + 1, en, a);
        goto loopend;
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
        // If we're at U2, don't allow leaving as nothing.
        auto right_unpaired = dp[piv + 1][en][DP_U];
        if (a != DP_U2) right_unpaired = std::min(right_unpaired, ZERO_E);

        // Check a == U_RC:
        // (   )<.( ** ). > Right coax backward
        if (a == DP_U_RC) {
          if (base11 + m->MismatchCoaxial(pl1b, pb, stb, st1b) + right_unpaired ==
              dp[st][en][DP_U_RC]) {
            // Ctds were already set from the recurrence that called this.
            q.emplace_back(st + 1, piv - 1, DP_P);
            if (right_unpaired != ZERO_E) q.emplace_back(piv + 1, en, DP_U);
            goto loopend;
          }
        }

        // DP_U_RC is only the above case.
        if (a == DP_U_RC) continue;

        // (   )<   > - U, U2, U_WC?, U_GU?
        if (base00 + right_unpaired == dp[st][en][a] && (a != DP_U_WC || IsWcPair(stb, pb)) &&
            (a != DP_U_GU || IsGuPair(stb, pb))) {
          // If U_WC, or U_GU, we were involved in some sort of coaxial stack previously, and were
          // already set.
          if (a != DP_U_WC && a != DP_U_GU) res.ctd[st] = CTD_UNUSED;
          q.emplace_back(st, piv, DP_P);
          if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace_back(piv + 1, en, DP_U);
          goto loopend;
        }

        // The rest of the cases are for U and U2.
        if (a != DP_U && a != DP_U2) continue;

        // (   )3<   > 3' - U, U2
        if (base01 + m->dangle3[pl1b][pb][stb] + right_unpaired == dp[st][en][a]) {
          res.ctd[st] = CTD_3_DANGLE;
          q.emplace_back(st, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace_back(piv + 1, en, DP_U);
          goto loopend;
        }
        // 5(   )<   > 5' - U, U2
        if (base10 + m->dangle5[pb][stb][st1b] + right_unpaired == dp[st][en][a]) {
          res.ctd[st + 1] = CTD_5_DANGLE;
          q.emplace_back(st + 1, piv, DP_P);
          if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace_back(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch - U, U2
        if (base11 + m->terminal[pl1b][pb][stb][st1b] + right_unpaired == dp[st][en][a]) {
          res.ctd[st + 1] = CTD_MISMATCH;
          q.emplace_back(st + 1, piv - 1, DP_P);
          if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace_back(piv + 1, en, DP_U);
          goto loopend;
        }
        // .(   ).<(   ) > Left coax - U, U2
        auto val = base11 + m->MismatchCoaxial(pl1b, pb, stb, st1b);
        if (val + dp[piv + 1][en][DP_U_WC] == dp[st][en][a]) {
          res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
          res.ctd[piv + 1] = CTD_LCOAX_WITH_PREV;
          q.emplace_back(st + 1, piv - 1, DP_P);
          q.emplace_back(piv + 1, en, DP_U_WC);
          goto loopend;
        }
        if (val + dp[piv + 1][en][DP_U_GU] == dp[st][en][a]) {
          res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
          res.ctd[piv + 1] = CTD_LCOAX_WITH_PREV;
          q.emplace_back(st + 1, piv - 1, DP_P);
          q.emplace_back(piv + 1, en, DP_U_GU);
          goto loopend;
        }

        // (   )<.(   ). > Right coax forward - U, U2
        if (base00 + dp[piv + 1][en][DP_U_RC] == dp[st][en][a]) {
          res.ctd[st] = CTD_RC_WITH_NEXT;
          res.ctd[piv + 2] = CTD_RC_WITH_PREV;
          q.emplace_back(st, piv, DP_P);
          q.emplace_back(piv + 1, en, DP_U_RC);
          goto loopend;
        }

        // (   )(<   ) > Flush coax - U, U2
        if (base01 + m->stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC] == dp[st][en][a]) {
          res.ctd[st] = CTD_FCOAX_WITH_NEXT;
          res.ctd[piv] = CTD_FCOAX_WITH_PREV;
          q.emplace_back(st, piv - 1, DP_P);
          q.emplace_back(piv, en, DP_U_WC);
          goto loopend;
        }
        if ((IsGu(pb)) &&
            base01 + m->stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU] == dp[st][en][a]) {
          res.ctd[st] = CTD_FCOAX_WITH_NEXT;
          res.ctd[piv] = CTD_FCOAX_WITH_PREV;
          q.emplace_back(st, piv - 1, DP_P);
          q.emplace_back(piv, en, DP_U_GU);
          goto loopend;
        }
      }
    }

  loopend: {};
  }

  return res;
}

}  // namespace mrna::md::base::opt
