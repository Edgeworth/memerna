// Copyright 2026 Eliot Courtney.
#include "backends/canon/trace/trace.h"

#include <algorithm>
#include <string>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "api/trace/trace_cfg.h"
#include "backends/canon/dp.h"
#include "backends/canon/energy/model.h"
#include "backends/common/static.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/error.h"
#include "util/log.h"
#include "util/util.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
bool Trace<CFG, PseudofreeCfgT>::IsSupported(
    const PseudofreeCfgT& /*pf*/, const trace::TraceCfg& tcfg, std::string* reason) {
  if (tcfg.random) {
    if (reason) *reason = "random traceback is not supported";
    return false;
  }
  return true;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
TraceResult Trace<CFG, PseudofreeCfgT>::Run(const Primary& r, const Model::Ptr& m,
    const DpState& state, const PseudofreeCfgT& pf, const trace::TraceCfg& tcfg) {
  const int N = r.size();
  tcfg.Verify();
  std::string reason;
  verify(IsSupported(pf, tcfg, &reason), "{} does not support the given configuration: {}",
      funcname(), reason);

  logdebug("canon {} with {}, {}, {}", funcname(), CFG, tcfg, pf);

  const auto& [dp, ext] = state;
  TraceResult res((Secondary(N)), Ctds(N));
  std::vector<DpIndex> q;
  q.emplace_back(0, INVALID_INDEX, EXT);
  while (!q.empty()) {
    const auto back = q.back();
    const int st = back.st;
    int en = back.en;
    const DpArrayId a = back.a;
    q.pop_back();

    if (back.IsExternal()) {
      if (a == EXT && st + 1 < N && ext[st + 1][EXT] + pf.Unpaired(st) == ext[st][EXT]) {
        q.emplace_back(st + 1, INVALID_INDEX, EXT);
        goto loopend;
      }
      if (a != EXT) continue;
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

        // (   )<   >
        auto val = base00 + ext[en + 1][EXT];
        Ctd val_ctd = CTD_UNUSED;
        if constexpr (CFG.UseD2()) {
          if (st != 0 && en != N - 1) {
            // (   )<   > Terminal mismatch
            val += m->terminal[enb][r[en + 1]][r[st - 1]][stb];
            val_ctd = CTD_MISMATCH;
          } else if (en != N - 1) {
            // (   )<3   > 3'
            val += m->dangle3[enb][r[en + 1]][stb];
            val_ctd = CTD_3_DANGLE;
          } else if (st != 0) {
            // 5(   )<   > 5'
            val += m->dangle5[enb][r[st - 1]][stb];
            val_ctd = CTD_5_DANGLE;
          }
        }

        if (val == ext[st][EXT]) {
          res.ctd[st] = val_ctd;
          q.emplace_back(st, en, DP_P);
          q.emplace_back(en + 1, INVALID_INDEX, EXT);
          goto loopend;
        }

        if constexpr (CFG.UseDangleMismatch()) {
          // (   )3<   > 3'
          if (base01 + m->dangle3[en1b][enb][stb] + pf.Unpaired(en) + ext[en + 1][EXT] ==
              ext[st][EXT]) {
            res.ctd[st] = CTD_3_DANGLE;
            q.emplace_back(st, en - 1, DP_P);
            q.emplace_back(en + 1, INVALID_INDEX, EXT);
            goto loopend;
          }
          // 5(   )<   > 5'
          if (base10 + m->dangle5[enb][stb][st1b] + pf.Unpaired(st) + ext[en + 1][EXT] ==
              ext[st][EXT]) {
            res.ctd[st + 1] = CTD_5_DANGLE;
            q.emplace_back(st + 1, en, DP_P);
            q.emplace_back(en + 1, INVALID_INDEX, EXT);
            goto loopend;
          }
          // .(   ).<   > Terminal mismatch
          if (base11 + m->terminal[en1b][enb][stb][st1b] + pf.Unpaired(st) + pf.Unpaired(en) +
                  ext[en + 1][EXT] ==
              ext[st][EXT]) {
            res.ctd[st + 1] = CTD_MISMATCH;
            q.emplace_back(st + 1, en - 1, DP_P);
            q.emplace_back(en + 1, INVALID_INDEX, EXT);
            goto loopend;
          }
        }

        if constexpr (CFG.UseCoaxialStacking()) {
          // .(   ).<(   ) > Left coax
          for (int piv = en + HAIRPIN_MIN_SZ + 2; piv < N; ++piv) {
            const Base prb = r[en + 1];
            const Base pb = r[piv];
            if (base11 + dp[en + 1][piv][DP_P] + m->AuGuPenalty(prb, pb) +
                    m->MismatchCoaxial(en1b, enb, stb, st1b) + pf.Unpaired(st) + pf.Unpaired(en) +
                    ext[piv + 1][EXT] ==
                ext[st][EXT]) {
              res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
              res.ctd[en + 1] = CTD_LCOAX_WITH_PREV;
              q.emplace_back(st + 1, en - 1, DP_P);
              q.emplace_back(en + 1, piv, DP_P);
              q.emplace_back(piv + 1, INVALID_INDEX, EXT);
              goto loopend;
            }
          }

          // (   )<.(   ). > Right coax forward
          for (int piv = en + HAIRPIN_MIN_SZ + 4; piv < N; ++piv) {
            const Base prb = r[en + 1];
            const Base pr1b = r[en + 2];
            const Base pl1b = r[piv - 1];
            const Base pb = r[piv];
            if (base00 + dp[en + 2][piv - 1][DP_P] + m->AuGuPenalty(pr1b, pl1b) +
                    m->MismatchCoaxial(pl1b, pb, prb, pr1b) + pf.Unpaired(en + 1) +
                    pf.Unpaired(piv) + ext[piv + 1][EXT] ==
                ext[st][EXT]) {
              res.ctd[st] = CTD_RC_WITH_NEXT;
              res.ctd[en + 2] = CTD_RC_WITH_PREV;
              q.emplace_back(st, en, DP_P);
              q.emplace_back(en + 2, piv - 1, DP_P);
              q.emplace_back(piv + 1, INVALID_INDEX, EXT);
              goto loopend;
            }
          }

          // (   )(<   ) > Flush coax
          for (int piv = en + HAIRPIN_MIN_SZ + 1; piv < N; ++piv) {
            const Base pb = r[piv];
            if (base01 + dp[en][piv][DP_P] + m->AuGuPenalty(enb, pb) +
                    m->stack[en1b][enb][pb][stb] + ext[piv + 1][EXT] ==
                ext[st][EXT]) {
              res.ctd[st] = CTD_FCOAX_WITH_NEXT;
              res.ctd[en] = CTD_FCOAX_WITH_PREV;
              q.emplace_back(st, en - 1, DP_P);
              q.emplace_back(en, piv, DP_P);
              q.emplace_back(piv + 1, INVALID_INDEX, EXT);
              goto loopend;
            }
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
        res.s[st] = As<Index>(en);
        res.s[en] = As<Index>(st);

        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (dp[ist][ien][DP_P] < CAP_E) {
              const auto val = m->TwoLoop<CFG>(r, pf, st, en, ist, ien) + dp[ist][ien][DP_P];
              if (val == dp[st][en][DP_P]) {
                q.emplace_back(ist, ien, DP_P);
                goto loopend;
              }
            }
          }
        }

        if (m->Hairpin<CFG>(r, pf, st, en) == dp[st][en][DP_P]) goto loopend;

        const auto base_branch_cost =
            m->AuGuPenalty(stb, enb) + pf.Paired(st, en) + m->multiloop_a + m->multiloop_b;

        // (<   ><   >)
        auto val = base_branch_cost + dp[st + 1][en - 1][DP_U2];
        Ctd val_ctd = CTD_UNUSED;
        if constexpr (CFG.UseD2()) {
          // D2 can overlap terminal mismatches with anything.
          // (<   ><   >) Terminal mismatch
          val += m->terminal[stb][st1b][en1b][enb];
          val_ctd = CTD_MISMATCH;
        }
        if (val == dp[st][en][DP_P]) {
          res.ctd[en] = val_ctd;
          q.emplace_back(st + 1, en - 1, DP_U2);
          goto loopend;
        }

        if constexpr (CFG.UseDangleMismatch()) {
          // (3<   ><   >) 3'
          if (base_branch_cost + dp[st + 2][en - 1][DP_U2] + m->dangle3[stb][st1b][enb] +
                  pf.Unpaired(st + 1) + m->multiloop_c ==
              dp[st][en][DP_P]) {
            res.ctd[en] = CTD_3_DANGLE;
            q.emplace_back(st + 2, en - 1, DP_U2);
            goto loopend;
          }
          // (<   ><   >5) 5'
          if (base_branch_cost + dp[st + 1][en - 2][DP_U2] + m->dangle5[stb][en1b][enb] +
                  pf.Unpaired(en - 1) + m->multiloop_c ==
              dp[st][en][DP_P]) {
            res.ctd[en] = CTD_5_DANGLE;
            q.emplace_back(st + 1, en - 2, DP_U2);
            goto loopend;
          }
          // (.<   ><   >.) Terminal mismatch
          if (base_branch_cost + dp[st + 2][en - 2][DP_U2] + m->terminal[stb][st1b][en1b][enb] +
                  pf.Unpaired(st + 1) + pf.Unpaired(en - 1) + 2 * m->multiloop_c ==
              dp[st][en][DP_P]) {
            res.ctd[en] = CTD_MISMATCH;
            q.emplace_back(st + 2, en - 2, DP_U2);
            goto loopend;
          }
        }

        if constexpr (CFG.UseCoaxialStacking()) {
          const auto outer_coax = m->MismatchCoaxial(stb, st1b, en1b, enb) + pf.Unpaired(st + 1) +
              pf.Unpaired(en - 1) + 2 * m->multiloop_c;
          for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
            const Base pl1b = r[piv - 1];
            const Base plb = r[piv];
            const Base prb = r[piv + 1];
            const Base pr1b = r[piv + 2];

            // (.(   )   .) Left outer coax - P
            if (base_branch_cost + dp[st + 2][piv][DP_P] + m->AuGuPenalty(st2b, plb) +
                    m->multiloop_b + dp[piv + 1][en - 2][DP_U] + outer_coax ==
                dp[st][en][DP_P]) {
              res.ctd[en] = CTD_LCOAX_WITH_NEXT;
              res.ctd[st + 2] = CTD_LCOAX_WITH_PREV;
              q.emplace_back(st + 2, piv, DP_P);
              q.emplace_back(piv + 1, en - 2, DP_U);
              goto loopend;
            }
            // (.   (   ).) Right outer coax
            if (base_branch_cost + dp[st + 2][piv][DP_U] + m->AuGuPenalty(prb, en2b) +
                    m->multiloop_b + dp[piv + 1][en - 2][DP_P] + outer_coax ==
                dp[st][en][DP_P]) {
              res.ctd[en] = CTD_RC_WITH_PREV;
              res.ctd[piv + 1] = CTD_RC_WITH_NEXT;
              q.emplace_back(st + 2, piv, DP_U);
              q.emplace_back(piv + 1, en - 2, DP_P);
              goto loopend;
            }

            // (.(   ).   ) Left inner coax
            if (base_branch_cost + dp[st + 2][piv - 1][DP_P] + m->AuGuPenalty(st2b, pl1b) +
                    m->multiloop_b + dp[piv + 1][en - 1][DP_U] +
                    m->MismatchCoaxial(pl1b, plb, st1b, st2b) + pf.Unpaired(st + 1) +
                    pf.Unpaired(piv) + 2 * m->multiloop_c ==
                dp[st][en][DP_P]) {
              res.ctd[en] = CTD_RC_WITH_NEXT;
              res.ctd[st + 2] = CTD_RC_WITH_PREV;
              q.emplace_back(st + 2, piv - 1, DP_P);
              q.emplace_back(piv + 1, en - 1, DP_U);
              goto loopend;
            }
            // (   .(   ).) Right inner coax
            if (base_branch_cost + dp[st + 1][piv][DP_U] + m->AuGuPenalty(pr1b, en2b) +
                    m->multiloop_b + dp[piv + 2][en - 2][DP_P] +
                    m->MismatchCoaxial(en2b, en1b, prb, pr1b) + pf.Unpaired(piv + 1) +
                    pf.Unpaired(en - 1) + 2 * m->multiloop_c ==
                dp[st][en][DP_P]) {
              res.ctd[en] = CTD_LCOAX_WITH_PREV;
              res.ctd[piv + 2] = CTD_LCOAX_WITH_NEXT;
              q.emplace_back(st + 1, piv, DP_U);
              q.emplace_back(piv + 2, en - 2, DP_P);
              goto loopend;
            }

            // ((   )   ) Left flush coax
            if (base_branch_cost + dp[st + 1][piv][DP_P] + m->AuGuPenalty(st1b, plb) +
                    m->multiloop_b + dp[piv + 1][en - 1][DP_U] + m->stack[stb][st1b][plb][enb] ==
                dp[st][en][DP_P]) {
              res.ctd[en] = CTD_FCOAX_WITH_NEXT;
              res.ctd[st + 1] = CTD_FCOAX_WITH_PREV;
              q.emplace_back(st + 1, piv, DP_P);
              q.emplace_back(piv + 1, en - 1, DP_U);
              goto loopend;
            }
            // (   (   )) Right flush coax
            if (base_branch_cost + dp[st + 1][piv][DP_U] + m->AuGuPenalty(prb, en1b) +
                    m->multiloop_b + dp[piv + 1][en - 1][DP_P] + m->stack[stb][prb][en1b][enb] ==
                dp[st][en][DP_P]) {
              res.ctd[en] = CTD_FCOAX_WITH_PREV;
              res.ctd[piv + 1] = CTD_FCOAX_WITH_NEXT;
              q.emplace_back(st + 1, piv, DP_U);
              q.emplace_back(piv + 1, en - 1, DP_P);
              goto loopend;
            }
          }
        }

        continue;
      }

      if (a == DP_CX) {
        if constexpr (CFG.UseCoaxialStacking()) {
          // (   )(<   ) Flush coax - CX
          for (int piv = st + HAIRPIN_MIN_SZ + 2; piv <= en - HAIRPIN_MIN_SZ - 1; ++piv) {
            const Base pl1b = r[piv - 1];
            const Base pb = r[piv];
            if (dp[st][piv - 1][DP_P] + m->AuGuPenalty(stb, pl1b) + m->multiloop_b +
                    dp[piv][en][DP_P] + m->AuGuPenalty(pb, enb) + m->multiloop_b +
                    m->stack[pl1b][pb][enb][stb] ==
                dp[st][en][DP_CX]) {
              res.ctd[st] = CTD_FCOAX_WITH_NEXT;
              res.ctd[piv] = CTD_FCOAX_WITH_PREV;
              q.emplace_back(st, piv - 1, DP_P);
              q.emplace_back(piv, en, DP_P);
              goto loopend;
            }
          }
          // .(   ).<(   ) Left coax - CX
          for (int piv = st + HAIRPIN_MIN_SZ + 3; piv <= en - HAIRPIN_MIN_SZ - 2; ++piv) {
            const Base pl1b = r[piv - 1];
            const Base pb = r[piv];
            const Base prb = r[piv + 1];
            if (dp[st + 1][piv - 1][DP_P] + m->AuGuPenalty(st1b, pl1b) + m->multiloop_b +
                    dp[piv + 1][en][DP_P] + m->AuGuPenalty(prb, enb) + m->multiloop_b +
                    m->MismatchCoaxial(pl1b, pb, stb, st1b) + pf.Unpaired(st) + pf.Unpaired(piv) +
                    2 * m->multiloop_c ==
                dp[st][en][DP_CX]) {
              res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
              res.ctd[piv + 1] = CTD_LCOAX_WITH_PREV;
              q.emplace_back(st + 1, piv - 1, DP_P);
              q.emplace_back(piv + 1, en, DP_P);
              goto loopend;
            }
          }
          // (   )<.(   ). Right coax - CX
          for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en - HAIRPIN_MIN_SZ - 4; ++piv) {
            const Base pb = r[piv];
            const Base prb = r[piv + 1];
            const Base pr1b = r[piv + 2];
            if (dp[st][piv][DP_P] + m->AuGuPenalty(stb, pb) + m->multiloop_b +
                    dp[piv + 2][en - 1][DP_P] + m->AuGuPenalty(pr1b, en1b) + m->multiloop_b +
                    m->MismatchCoaxial(en1b, enb, prb, pr1b) + pf.Unpaired(piv + 1) +
                    pf.Unpaired(en) + 2 * m->multiloop_c ==
                dp[st][en][DP_CX]) {
              res.ctd[st] = CTD_RC_WITH_NEXT;
              res.ctd[piv + 2] = CTD_RC_WITH_PREV;
              q.emplace_back(st, piv, DP_P);
              q.emplace_back(piv + 2, en - 1, DP_P);
              goto loopend;
            }
          }
        }
        continue;
      }

      if (st + 1 < en && (a == DP_U || a == DP_U2) &&
          dp[st + 1][en][a] + pf.Unpaired(st) + m->multiloop_c == dp[st][en][a]) {
        q.emplace_back(st + 1, en, a);
        goto loopend;
      }

      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = r[piv];
        const auto pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on
        // the right.
        const auto base00 = dp[st][piv][DP_P] + m->AuGuPenalty(stb, pb) + m->multiloop_b;
        const auto base01 = dp[st][piv - 1][DP_P] + m->AuGuPenalty(stb, pl1b) + m->multiloop_b;
        const auto base10 = dp[st + 1][piv][DP_P] + m->AuGuPenalty(st1b, pb) + m->multiloop_b;
        const auto base11 = dp[st + 1][piv - 1][DP_P] + m->AuGuPenalty(st1b, pl1b) + m->multiloop_b;

        auto right_unpaired = dp[piv + 1][en][DP_U];
        bool is_right_filled = true;
        if (a != DP_U2 &&
            pf.UnpairedSum(piv + 1, en) + (en - piv) * m->multiloop_c < right_unpaired) {
          right_unpaired = pf.UnpairedSum(piv + 1, en) + (en - piv) * m->multiloop_c;
          is_right_filled = false;
        }

        // (   )<   > - U, U2
        auto val = base00 + right_unpaired;
        Ctd val_ctd = CTD_UNUSED;

        if constexpr (CFG.UseD2()) {
          // Note that D2 can overlap with anything.
          if (st != 0 && piv != N - 1) {
            // (   )<   > Terminal mismatch - U, U2
            val += m->terminal[pb][r[piv + 1]][r[st - 1]][stb];
            val_ctd = CTD_MISMATCH;
          } else if (piv != N - 1) {
            // (   )<3   > 3' - U, U2
            val += m->dangle3[pb][r[piv + 1]][stb];
            val_ctd = CTD_3_DANGLE;
          } else if (st != 0) {
            // 5(   )<   > 5' - U, U2
            val += m->dangle5[pb][r[st - 1]][stb];
            val_ctd = CTD_5_DANGLE;
          }
        }

        if (val == dp[st][en][a]) {
          res.ctd[st] = val_ctd;
          q.emplace_back(st, piv, DP_P);
          if (is_right_filled) q.emplace_back(piv + 1, en, DP_U);
          goto loopend;
        }

        if constexpr (CFG.UseDangleMismatch()) {
          // (   )3<   > 3' - U, U2
          if (base01 + m->dangle3[pl1b][pb][stb] + pf.Unpaired(piv) + m->multiloop_c +
                  right_unpaired ==
              dp[st][en][a]) {
            res.ctd[st] = CTD_3_DANGLE;
            q.emplace_back(st, piv - 1, DP_P);
            if (is_right_filled) q.emplace_back(piv + 1, en, DP_U);
            goto loopend;
          }
          // 5(   )<   > 5' - U, U2
          if (base10 + m->dangle5[pb][stb][st1b] + pf.Unpaired(st) + m->multiloop_c +
                  right_unpaired ==
              dp[st][en][a]) {
            res.ctd[st + 1] = CTD_5_DANGLE;
            q.emplace_back(st + 1, piv, DP_P);
            if (is_right_filled) q.emplace_back(piv + 1, en, DP_U);
            goto loopend;
          }
          // .(   ).<   > Terminal mismatch - U, U2
          if (base11 + m->terminal[pl1b][pb][stb][st1b] + pf.Unpaired(st) + pf.Unpaired(piv) +
                  2 * m->multiloop_c + right_unpaired ==
              dp[st][en][a]) {
            res.ctd[st + 1] = CTD_MISMATCH;
            q.emplace_back(st + 1, piv - 1, DP_P);
            if (is_right_filled) q.emplace_back(piv + 1, en, DP_U);
            goto loopend;
          }
        }
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        if (a == DP_U || a == DP_U2) {
          for (int piv = st + 2 * HAIRPIN_MIN_SZ + 3; piv <= en; ++piv) {
            // (   )(<   )<   > - CX
            auto right_unpaired = dp[piv + 1][en][DP_U];
            bool is_right_filled = true;
            if (pf.UnpairedSum(piv + 1, en) + (en - piv) * m->multiloop_c < right_unpaired) {
              right_unpaired = pf.UnpairedSum(piv + 1, en) + (en - piv) * m->multiloop_c;
              is_right_filled = false;
            }
            if (dp[st][piv][DP_CX] + right_unpaired == dp[st][en][a]) {
              q.emplace_back(st, piv, DP_CX);
              if (is_right_filled) q.emplace_back(piv + 1, en, DP_U);
              goto loopend;
            }
          }
        }
      }
    }

  loopend: {};
  }

  return res;
}

#define MRNA_INSTANTIATE(CFG, PF) template class Trace<CFG, PF>;
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::canon
