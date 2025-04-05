// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <utility>

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/precomp.h"
#include "backends/base/pfn/pfn.h"
#include "backends/common/base/dp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/pfn.h"
#include "model/primary.h"
#include "util/error.h"
#include "util/util.h"

namespace mrna::md::base {

PfnTables PfnDebug(const Primary& r, const Model::Ptr& initial_m, PfnState& state) {
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  // Force bulge states off.
  auto m = initial_m->Clone();
  auto cfg = m->cfg();
  cfg.bulge_states = false;
  m->SetEnergyCfg(cfg);

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false},  // Bulge states with partition function doesn't make sense.
      .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::D2,
          erg::EnergyCfg::Ctd::NONE},
  };
  support.VerifySupported(funcname(), m->cfg());
  m->pf.Verify(r);

  spdlog::debug("base {} with cfg {}", funcname(), m->cfg());

  const int N = static_cast<int>(r.size());
  const Precomp pc(Primary(r), m);
  state.dp = BoltzDpArray(r.size() + 1, 0);
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
        BoltzEnergy p{0};
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist)
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien)
            p += pc.TwoLoop(st, en, ist, ien).Boltz() * dp[ist][ien][PT_P];
        // Hairpin loops.
        p += m->Hairpin(r, st, en).Boltz();

        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const BoltzEnergy base_branch_cost =
            (pc.augubranch[stb][enb] + m->pf.Paired(st, en) + m->multiloop_hack_a).Boltz();

        // (<   ><   >)
        BoltzEnergy val = base_branch_cost * dp[st + 1][en - 1][PT_U2];
        if (m->cfg().UseD2()) {
          // D2 can overlap terminal mismatches with anything.
          // (<   ><   >) Terminal mismatch
          val *= m->terminal[stb][st1b][en1b][enb].Boltz();
        }
        p += val;

        if (m->cfg().UseDangleMismatch()) {
          // (3<   ><   >) 3'
          p += base_branch_cost * dp[st + 2][en - 1][PT_U2] *
              (m->dangle3[stb][st1b][enb] + m->pf.Unpaired(st + 1)).Boltz();
          // (<   ><   >5) 5'
          p += base_branch_cost * dp[st + 1][en - 2][PT_U2] *
              (m->dangle5[stb][en1b][enb] + m->pf.Unpaired(en - 1)).Boltz();
          // (.<   ><   >.) Terminal mismatch
          p += base_branch_cost * dp[st + 2][en - 2][PT_U2] *
              (m->terminal[stb][st1b][en1b][enb] + m->pf.Unpaired(st + 1) + m->pf.Unpaired(en - 1))
                  .Boltz();
        }

        if (m->cfg().UseCoaxialStacking()) {
          const auto outer_coax = m->MismatchCoaxial(stb, st1b, en1b, enb) +
              m->pf.Unpaired(st + 1) + m->pf.Unpaired(en - 1);
          for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
            // Paired coaxial stacking cases:
            const Base pl1b = r[piv - 1];
            const Base plb = r[piv];
            const Base prb = r[piv + 1];
            const Base pr1b = r[piv + 2];

            // (.(   )   .) Left outer coax - P
            p += base_branch_cost * dp[st + 2][piv][PT_P] * dp[piv + 1][en - 2][PT_U] *
                (pc.augubranch[st2b][plb] + outer_coax).Boltz();
            // (.   (   ).) Right outer coax
            p += base_branch_cost * dp[st + 2][piv][PT_U] * dp[piv + 1][en - 2][PT_P] *
                (pc.augubranch[prb][en2b] + outer_coax).Boltz();

            // (.(   ).   ) Left inner coax
            p += base_branch_cost * dp[st + 2][piv - 1][PT_P] * dp[piv + 1][en - 1][PT_U] *
                (pc.augubranch[st2b][pl1b] + m->MismatchCoaxial(pl1b, plb, st1b, st2b) +
                    m->pf.Unpaired(st + 1) + m->pf.Unpaired(piv))
                    .Boltz();
            // (   .(   ).) Right inner coax
            p += base_branch_cost * dp[st + 1][piv][PT_U] * dp[piv + 2][en - 2][PT_P] *
                (pc.augubranch[pr1b][en2b] + m->MismatchCoaxial(en2b, en1b, prb, pr1b) +
                    m->pf.Unpaired(piv + 1) + m->pf.Unpaired(en - 1))
                    .Boltz();

            // ((   )   ) Left flush coax
            p += base_branch_cost * dp[st + 1][piv][PT_P] * dp[piv + 1][en - 1][PT_U] *
                (pc.augubranch[st1b][plb] + m->stack[stb][st1b][plb][enb]).Boltz();
            // (   (   )) Right flush coax
            p += base_branch_cost * dp[st + 1][piv][PT_U] * dp[piv + 1][en - 1][PT_P] *
                (pc.augubranch[prb][en1b] + m->stack[stb][prb][en1b][enb]).Boltz();
          }
        }

        dp[st][en][PT_P] = p;
      }
      BoltzEnergy u{0};
      BoltzEnergy u2{0};
      BoltzEnergy rcoax{0};
      BoltzEnergy wc{0};
      BoltzEnergy gu{0};
      // Update unpaired.
      // Choose `st` to be unpaired.
      if (st + 1 < en) {
        u += dp[st + 1][en][PT_U] * m->pf.Unpaired(st).Boltz();
        u2 += dp[st + 1][en][PT_U2] * m->pf.Unpaired(st).Boltz();
      }

      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = r[piv];
        const auto pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const BoltzEnergy base00 = dp[st][piv][PT_P] * pc.augubranch[stb][pb].Boltz();
        const BoltzEnergy base01 = dp[st][piv - 1][PT_P] * pc.augubranch[stb][pl1b].Boltz();
        const BoltzEnergy base10 = dp[st + 1][piv][PT_P] * pc.augubranch[st1b][pb].Boltz();
        const BoltzEnergy base11 = dp[st + 1][piv - 1][PT_P] * pc.augubranch[st1b][pl1b].Boltz();

        const BoltzEnergy right_paired = dp[piv + 1][en][PT_U];
        const BoltzEnergy right_unpaired = m->pf.UnpairedCum(piv + 1, en).Boltz();

        // (   )<   > - U, U_WC?, U_GU?
        BoltzEnergy u2_val = base00 * right_paired;
        BoltzEnergy val = base00 * right_unpaired + base00 * right_paired;

        if (m->cfg().UseD2()) {
          // Note that D2 can overlap with anything.
          if (st != 0 && piv != N - 1) {
            // (   )<   > Terminal mismatch - U
            val *= m->terminal[pb][r[piv + 1]][r[st - 1]][stb].Boltz();
            u2_val *= m->terminal[pb][r[piv + 1]][r[st - 1]][stb].Boltz();
          } else if (piv != N - 1) {
            // (   )<3   > 3' - U
            val *= m->dangle3[pb][r[piv + 1]][stb].Boltz();
            u2_val *= m->dangle3[pb][r[piv + 1]][stb].Boltz();
          } else if (st != 0) {
            // 5(   )<   > 5' - U
            val *= m->dangle5[pb][r[st - 1]][stb].Boltz();
            u2_val *= m->dangle5[pb][r[st - 1]][stb].Boltz();
          }
        }

        u2 += u2_val;
        u += val;
        if (IsGuPair(stb, pb))
          gu += val;
        else
          wc += val;

        if (m->cfg().UseDangleMismatch()) {
          // (   )3<   > 3' - U
          val = base01 * (m->dangle3[pl1b][pb][stb] + m->pf.Unpaired(piv)).Boltz();
          u += val * right_unpaired + val * right_paired;
          u2 += val * right_paired;

          // 5(   )<   > 5' - U
          val = base10 * (m->dangle5[pb][stb][st1b] + m->pf.Unpaired(st)).Boltz();
          u += val * right_unpaired + val * right_paired;
          u2 += val * right_paired;

          // .(   ).<   > Terminal mismatch - U
          val = base11 *
              (m->terminal[pl1b][pb][stb][st1b] + m->pf.Unpaired(st) + m->pf.Unpaired(piv)).Boltz();
          u += val * right_unpaired + val * right_paired;
          u2 += val * right_paired;
        }

        if (m->cfg().UseCoaxialStacking()) {
          // .(   ).<(   ) > Left coax - U
          val = base11 *
              (m->MismatchCoaxial(pl1b, pb, stb, st1b) + m->pf.Unpaired(st) + m->pf.Unpaired(piv))
                  .Boltz() *
              (dp[piv + 1][en][PT_U_WC] + dp[piv + 1][en][PT_U_GU]);
          u += val;
          u2 += val;

          // (   )<.(   ). > Right coax forward and backward
          val = base00 * dp[piv + 1][en][PT_U_RC];
          u += val;
          u2 += val;

          val = base11 *
              (m->MismatchCoaxial(pl1b, pb, stb, st1b) + m->pf.Unpaired(st) + m->pf.Unpaired(piv))
                  .Boltz();
          rcoax += val * right_unpaired + val * right_paired;

          // (   )(<   ) > Flush coax - U
          val = base01 * m->stack[pl1b][pb][WcPair(pb)][stb].Boltz() * dp[piv][en][PT_U_WC];
          u += val;
          u2 += val;
          if (IsGu(pb)) {
            val = base01 * m->stack[pl1b][pb][GuPair(pb)][stb].Boltz() * dp[piv][en][PT_U_GU];
            u += val;
            u2 += val;
          }
        }
      }
      dp[st][en][PT_U] = u;
      dp[st][en][PT_U2] = u2;
      dp[st][en][PT_U_WC] = wc;
      dp[st][en][PT_U_GU] = gu;
      dp[st][en][PT_U_RC] = rcoax;
    }
  }

  // Compute the exterior tables.
  PfnExterior(r, *m, state);
  const auto& ext = state.ext;

  // Fill the left triangle.
  // The meaning of the tables changes here:
  // U, U2, WC, GU, RC: any table index with en < st must have a loop enclosing (en, st)
  for (int st = N - 1; st >= 0; --st) {
    for (int en = 0; en < st; ++en) {
      //        ..)...(..
      // rspace  en   st  lspace
      const int lspace = N - st - 1;
      const int rspace = en;
      const Base stb = r[st];
      const Base st1b = lspace ? r[st + 1] : Base(-1);
      const Base st2b = lspace > 1 ? r[st + 2] : Base(-1);
      const Base enb = r[en];
      const Base en1b = rspace ? r[en - 1] : Base(-1);
      const Base en2b = rspace > 1 ? r[en - 2] : Base(-1);

      if (m->CanPair(r, en, st)) {
        BoltzEnergy p = ZERO_B;
        const int ost_max = std::min(st + TWOLOOP_MAX_SZ + 2, N);
        for (int ost = st + 1; ost < ost_max; ++ost) {
          const int oen_min = std::max(en - TWOLOOP_MAX_SZ - 1 + (ost - st - 1), 0);
          for (int oen = en - 1; oen >= oen_min; --oen) {
            p += pc.TwoLoop(oen, ost, en, st).Boltz() * dp[ost][oen][PT_P];
          }
        }
        const BoltzEnergy base_branch_cost =
            (pc.augubranch[stb][enb] + m->multiloop_hack_a).Boltz();
        const Energy outer_coax = lspace && rspace ? m->MismatchCoaxial(stb, st1b, en1b, enb) +
                m->pf.Unpaired(st + 1) + m->pf.Unpaired(en - 1)
                                                   : MAX_E;
        // Try being an exterior loop - coax cases handled in the loop after this.
        {
          const BoltzEnergy augu = m->AuGuPenalty(enb, stb).Boltz();
          const BoltzEnergy rext = ext[st + 1][PTEXT_R];
          const BoltzEnergy r1ext = lspace > 1 ? ext[st + 2][PTEXT_R] : BoltzEnergy{1};
          const BoltzEnergy lext = rspace ? ext[en - 1][PTEXT_L] : BoltzEnergy{1};
          const BoltzEnergy l1ext = rspace > 1 ? ext[en - 2][PTEXT_L] : BoltzEnergy{1};

          BoltzEnergy d2_val = ONE_B;
          if (m->cfg().UseD2()) {
            if (st != N - 1 && en != 0) {
              // |<   >)   (<   >| Terminal mismatch
              d2_val *= m->terminal[stb][st1b][en1b][enb].Boltz();
            } else if (st != N - 1) {
              // |<   >)   (<3   >| 3'
              d2_val *= m->dangle3[stb][st1b][enb].Boltz();
            } else if (en != 0) {
              // |<   5>)   (<   >| 5'
              d2_val *= m->dangle5[stb][en1b][enb].Boltz();
            }
          }
          // |<   >)   (<   >| - Exterior loop
          p += augu * lext * rext * d2_val;

          // |   >)   (<   | - Enclosing loop
          if (lspace && rspace) p += base_branch_cost * dp[st + 1][en - 1][PT_U2] * d2_val;

          if (m->cfg().UseDangleMismatch()) {
            if (lspace) {
              // |<   >(   )3<   >| 3' - Exterior loop
              // lspace > 0
              p += augu * lext * r1ext *
                  (m->dangle3[stb][st1b][enb] + m->pf.Unpaired(st + 1)).Boltz();
              // |  >5)   (<   | 5' - Enclosing loop
              if (rspace > 1)
                p += base_branch_cost * dp[st + 1][en - 2][PT_U2] *
                    (m->dangle5[stb][en1b][enb] + m->pf.Unpaired(en - 1)).Boltz();
            }
            if (rspace) {
              // |<   >5(   )<   >| 5' - Exterior loop
              // rspace > 0
              p += augu * l1ext * rext *
                  (m->dangle5[stb][en1b][enb] + m->pf.Unpaired(en - 1)).Boltz();
              // |   >)   (3<  | 3' - Enclosing loop
              if (lspace > 1)
                p += base_branch_cost * dp[st + 2][en - 1][PT_U2] *
                    (m->dangle3[stb][st1b][enb] + m->pf.Unpaired(st + 1)).Boltz();
            }
            // |<   >m(   )m<   >| Terminal mismatch - Exterior loop
            if (lspace && rspace)
              p += augu * l1ext * r1ext *
                  (m->terminal[stb][st1b][en1b][enb] + m->pf.Unpaired(st + 1) +
                      m->pf.Unpaired(en - 1))
                      .Boltz();

            // |  >m)   (m<  | Terminal mismatch - Enclosing loop
            if (lspace > 1 && rspace > 1)
              p += base_branch_cost * dp[st + 2][en - 2][PT_U2] *
                  (m->terminal[stb][st1b][en1b][enb] + m->pf.Unpaired(st + 1) +
                      m->pf.Unpaired(en - 1))
                      .Boltz();
          }

          if (m->cfg().UseCoaxialStacking()) {
            const int limit = en + N;
            for (int tpiv = st; tpiv <= limit; ++tpiv) {
              const int pl = FastMod(tpiv - 1, N);
              const int piv = FastMod(tpiv, N);
              const int pr = FastMod(tpiv + 1, N);
              const Base pl1b = r[pl];
              const Base plb = r[piv];
              const Base prb = r[pr];
              const bool left_formable = tpiv < N && tpiv - st - 2 >= HAIRPIN_MIN_SZ;
              const bool right_formable = tpiv >= N && en - piv - 2 >= HAIRPIN_MIN_SZ;

              if (left_formable) {
                const BoltzEnergy rp1ext = ext[piv + 1][PTEXT_R];
                // |<   >)   (.(   ).<   >| Exterior loop - Left inner coax
                // lspace > 1 && not enclosed
                p += PairedWithPf(m, dp, st + 2, pl) * lext * rp1ext *
                    (m->AuGuPenalty(stb, enb) + m->AuGuPenalty(st2b, pl1b) +
                        m->MismatchCoaxial(pl1b, plb, st1b, st2b) + m->pf.Unpaired(st + 1) +
                        m->pf.Unpaired(piv))
                        .Boltz();
                // |<   >)   ((   )<   >| Exterior loop - Left flush coax
                // lspace > 0 && not enclosed
                p += PairedWithPf(m, dp, st + 1, piv) * lext * rp1ext *
                    (m->AuGuPenalty(stb, enb) + m->AuGuPenalty(st1b, plb) +
                        m->stack[stb][st1b][plb][enb])
                        .Boltz();

                if (rspace) {
                  // |<   >.)   (.(   )<   >| Exterior loop - Left outer coax
                  // lspace > 0 && rspace > 0 && not enclosed
                  p += PairedWithPf(m, dp, st + 2, piv) * l1ext * rp1ext *
                      (m->AuGuPenalty(stb, enb) + m->AuGuPenalty(st2b, plb) + outer_coax).Boltz();
                }
              }

              if (right_formable) {
                const BoltzEnergy lpext = piv > 0 ? ext[piv - 1][PTEXT_L] : BoltzEnergy{1};
                // |<   >.(   ).)   (<   >| Exterior loop - Right inner coax
                // rspace > 1 && not enclosed
                p += PairedWithPf(m, dp, pr, en - 2) * lpext * rext *
                    (m->AuGuPenalty(stb, enb) + m->AuGuPenalty(prb, en2b) +
                        m->MismatchCoaxial(en2b, en1b, plb, prb) + m->pf.Unpaired(piv) +
                        m->pf.Unpaired(en - 1))
                        .Boltz();
                // |<   >(   ))   (<   >| Exterior loop - Right flush coax
                // rspace > 0 && not enclosed
                p += PairedWithPf(m, dp, piv, en - 1) * lpext * rext *
                    (m->AuGuPenalty(stb, enb) + m->AuGuPenalty(plb, en1b) +
                        m->stack[en1b][enb][stb][plb])
                        .Boltz();

                if (lspace) {
                  // |<   >(   ).)   (.<   >| Exterior loop - Right outer coax
                  // lspace > 0 && rspace > 1 && not enclosed
                  p += PairedWithPf(m, dp, piv, en - 2) * lpext * r1ext *
                      (m->AuGuPenalty(stb, enb) + m->AuGuPenalty(plb, en2b) + outer_coax).Boltz();
                }
              }
            }
          }
        }

        // Enclosing loop cases.
        // Can start at st + 2 because we need to form an enclosing loop.
        // At worst the enclosing loop has to start at st + 1.
        if (m->cfg().UseCoaxialStacking()) {
          const int limit = en + N;
          for (int tpiv = st + 2; tpiv < limit; ++tpiv) {
            const int pl = FastMod(tpiv - 1, N);
            const int piv = FastMod(tpiv, N);
            const int pr = FastMod(tpiv + 1, N);
            const int pr1 = FastMod(tpiv + 2, N);
            // Left block is: [st, piv], Right block is: [piv + 1, en].
            const Base pl1b = r[pl];
            const Base plb = r[piv];
            const Base prb = r[pr];
            const Base pr1b = r[pr1];
            // When neither the left block nor the right block straddles the border
            // we don't get enclosing loops sometimes, so check against N - 1.
            const bool straddling = tpiv != (N - 1);
            // Left loop formable if straddling, and is big enough or crosses over.
            const bool left_formable = straddling && (piv - st - 2 >= HAIRPIN_MIN_SZ || tpiv >= N);
            const bool right_formable =
                straddling && (en - piv - 3 >= HAIRPIN_MIN_SZ || tpiv < N - 1);
            const bool left_dot_formable = left_formable && tpiv != N;  // Can't split a dot.
            const bool right_dot_formable = right_formable && tpiv != N - 2;  // Can't split a dot.

            if (lspace > 1 && rspace > 1) {
              if (left_formable) {
                // |  >.)   (.(   )<  | Enclosing loop - Left outer coax
                // lspace > 1 && rspace > 1 && enclosed
                p += base_branch_cost * PairedWithPf(m, dp, st + 2, piv) * dp[pr][en - 2][PT_U] *
                    (pc.augubranch[st2b][plb] + outer_coax).Boltz();
              }

              if (right_formable) {
                // |  >(   ).)   (.<  | Enclosing loop - Right outer coax
                // lspace > 1 && rspace > 1 && enclosed
                p += base_branch_cost * dp[st + 2][piv][PT_U] * PairedWithPf(m, dp, pr, en - 2) *
                    (pc.augubranch[prb][en2b] + outer_coax).Boltz();
              }
            }

            if (lspace > 1 && rspace && left_dot_formable) {
              // |  >)   (.(   ).<  | Enclosing loop - Left inner coax
              // lspace > 1 && rspace > 0 && enclosed && no dot split
              p += base_branch_cost * PairedWithPf(m, dp, st + 2, pl) * dp[pr][en - 1][PT_U] *
                  (pc.augubranch[st2b][pl1b] + m->MismatchCoaxial(pl1b, plb, st1b, st2b) +
                      m->pf.Unpaired(st + 1) + m->pf.Unpaired(piv))
                      .Boltz();
            }

            if (lspace && rspace > 1 && right_dot_formable) {
              // |  >.(   ).)   (<  | Enclosing loop - Right inner coax
              // lspace > 0 && rspace > 1 && enclosed && no dot split
              p += base_branch_cost * dp[st + 1][piv][PT_U] * PairedWithPf(m, dp, pr1, en - 2) *
                  (pc.augubranch[pr1b][en2b] + m->MismatchCoaxial(en2b, en1b, prb, pr1b) +
                      m->pf.Unpaired(pr) + m->pf.Unpaired(en - 1))
                      .Boltz();
            }

            if (lspace && rspace) {
              if (left_formable) {
                // |  >)   ((   )<  | Enclosing loop - Left flush coax
                // lspace > 0 && rspace > 0 && enclosed
                p += base_branch_cost * PairedWithPf(m, dp, st + 1, piv) * dp[pr][en - 1][PT_U] *
                    (pc.augubranch[st1b][plb] + m->stack[stb][st1b][plb][enb]).Boltz();
              }

              if (right_formable) {
                // |  >(   ))   (<  | Enclosing loop - Right flush coax
                // lspace > 0 && rspace > 0 && enclosed
                p += base_branch_cost * dp[st + 1][piv][PT_U] * PairedWithPf(m, dp, pr, en - 1) *
                    (pc.augubranch[prb][en1b] + m->stack[en1b][enb][stb][prb]).Boltz();
              }
            }
          }
        }

        dp[st][en][PT_P] = p;
      }

      BoltzEnergy u{0};
      BoltzEnergy u2{0};
      BoltzEnergy rcoax{0};
      BoltzEnergy wc{0};
      BoltzEnergy gu{0};
      // Update unpaired.
      // Choose `st` to be unpaired, but only if we can maintain the constraint that we have
      // an enclosing loop formed.
      if (st + 1 < N) {
        u += dp[st + 1][en][PT_U] * m->pf.Unpaired(st).Boltz();
        u2 += dp[st + 1][en][PT_U2] * m->pf.Unpaired(st).Boltz();
      }

      for (int tpiv = st + 1; tpiv <= en + N; ++tpiv) {
        const int pl = FastMod(tpiv - 1, N);
        const int piv = FastMod(tpiv, N);
        const int pr = FastMod(tpiv + 1, N);
        const auto pb = r[piv];
        const auto pl1b = r[pl];
        const auto prb = r[pr];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const BoltzEnergy base00 = PairedWithPf(m, dp, st, piv) * pc.augubranch[stb][pb].Boltz();
        const BoltzEnergy base01 = PairedWithPf(m, dp, st, pl) * pc.augubranch[stb][pl1b].Boltz();

        const BoltzEnergy right_paired = dp[pr][en][PT_U];
        const BoltzEnergy right_unpaired = tpiv >= N ? m->pf.UnpairedCum(pr, en).Boltz() : ONE_B;

        // Must have an enclosing loop.
        const bool straddling = tpiv != N - 1;
        const bool dot_straddling = straddling && tpiv != N;
        BoltzEnergy val{0};

        if (straddling) {
          BoltzEnergy d2_val = ONE_B;
          if (m->cfg().UseD2()) {
            // |  >>   m<(   )<m  | Terminal mismatch
            // tpiv != N-1 && st != 0, so we can always do a terminal mismatch.
            d2_val *= m->terminal[pb][prb][r[st - 1]][stb].Boltz();
          }

          // |  >>   <(   )<  |
          val = base00 * right_paired * d2_val;
          u2 += val;
          u += val;
          if (IsGuPair(stb, pb))
            gu += val;
          else
            wc += val;

          // U must cross the boundary to have the rest of it be nothing.
          if (tpiv >= N) {
            val = base00 * right_unpaired * d2_val;
            u += val;
            if (IsGuPair(stb, pb))
              gu += val;
            else
              wc += val;
          }

          if (m->cfg().UseCoaxialStacking()) {
            // |  )  >>   <(   )<(  | Flush coax
            // straddling
            val = base00 * m->stack[pb][prb][WcPair(prb)][stb].Boltz() * dp[pr][en][PT_U_WC];
            u += val;
            u2 += val;
            if (IsGu(prb)) {
              val = base00 * m->stack[pb][prb][GuPair(prb)][stb].Boltz() * dp[pr][en][PT_U_GU];
              u += val;
              u2 += val;
            }
            // |  ).  >>   <(   )<.(   | Right coax forward
            // straddling
            val = base00 * dp[pr][en][PT_U_RC];
            u += val;
            u2 += val;
          }
        }

        if (m->cfg().UseDangleMismatch() && dot_straddling) {
          // |  >>   <(   )3<  | 3'
          val = base01 * (m->dangle3[pl1b][pb][stb] + m->pf.Unpaired(piv)).Boltz();
          if (tpiv >= N) u += val * right_unpaired;
          val *= right_paired;
          u += val;
          u2 += val;
        }

        if (lspace) {
          const BoltzEnergy base10 =
              PairedWithPf(m, dp, st + 1, piv) * pc.augubranch[st1b][pb].Boltz();
          const BoltzEnergy base11 =
              PairedWithPf(m, dp, st + 1, pl) * pc.augubranch[st1b][pl1b].Boltz();

          if (m->cfg().UseDangleMismatch() && straddling) {
            // |  >>   <5(   )<  | 5'
            val = base10 * (m->dangle5[pb][stb][st1b] + m->pf.Unpaired(st)).Boltz();
            if (tpiv >= N) u += val * right_unpaired;
            val *= right_paired;
            u += val;
            u2 += val;
          }

          if (m->cfg().UseDangleMismatch() && dot_straddling) {
            // |  >>   <m(   )m<  | Terminal mismatch
            // lspace > 0 && dot_straddling
            val = base11 *
                (m->terminal[pl1b][pb][stb][st1b] + m->pf.Unpaired(st) + m->pf.Unpaired(piv))
                    .Boltz();
            if (tpiv >= N) u += val * right_unpaired;
            val *= right_paired;
            u += val;
            u2 += val;
          }

          if (m->cfg().UseCoaxialStacking() && dot_straddling) {
            // |  )>>   <.(   ).<(  | Left coax
            // lspace > 0 && dot_straddling
            val = base11 *
                (m->MismatchCoaxial(pl1b, pb, stb, st1b) + m->pf.Unpaired(st) + m->pf.Unpaired(piv))
                    .Boltz();
            val = val * (dp[pr][en][PT_U_WC] + dp[pr][en][PT_U_GU]);
            u += val;
            u2 += val;
            // |  ).  >>   <(   )<.(   | Right coax backward
            val = base11 *
                (m->MismatchCoaxial(pl1b, pb, stb, st1b) + m->pf.Unpaired(st) + m->pf.Unpaired(piv))
                    .Boltz();
            if (tpiv >= N) rcoax += val * right_unpaired;
            rcoax += val * right_paired;
          }
        }
      }

      dp[st][en][PT_U] = u;
      dp[st][en][PT_U2] = u2;
      dp[st][en][PT_U_WC] = wc;
      dp[st][en][PT_U_GU] = gu;
      dp[st][en][PT_U_RC] = rcoax;
    }
  }

  BoltzSums p(N, 0);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) p[i][j] = dp[i][j][PT_P];

  return {std::move(p), ext[0][PTEXT_R]};
}

}  // namespace mrna::md::base
