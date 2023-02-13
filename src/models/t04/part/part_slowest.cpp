// Copyright 2016 Eliot Courtney.
#include <algorithm>
#include <memory>
#include <utility>

#include "api/energy/energy_cfg.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/part.h"
#include "model/primary.h"
#include "models/t04/energy/model.h"
#include "models/t04/energy/precomp.h"
#include "models/t04/part/part.h"
#include "util/util.h"

namespace mrna::md::t04 {

Part PartitionSlowest(const Primary& r, const Model::Ptr& initial_em, PartState& state) {
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  // Force bulge states off.
  auto em = initial_em->Clone();
  em->cfg.bulge_states = false;

  verify(em->cfg.lonely_pairs != erg::EnergyCfg::LonelyPairs::OFF,
      "fully disallowing lonely pairs is not supported in this energy model");
  verify(
      em->cfg.ctd == erg::EnergyCfg::Ctd::ALL, "only full CTDs are supported in this energy model");

  const int N = static_cast<int>(r.size());
  const Precomp pc(Primary(r), em);
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

      if (em->CanPair(r, st, en)) {
        BoltzEnergy p{0};
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist)
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien)
            p += pc.TwoLoop(st, en, ist, ien).Boltz() * dp[ist][ien][PT_P];
        // Hairpin loops.
        p += em->Hairpin(r, st, en).Boltz();
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const BoltzEnergy base_branch_cost =
            (pc.augubranch[stb][enb] + em->multiloop_hack_a).Boltz();

        // (<   ><   >)
        p += base_branch_cost * dp[st + 1][en - 1][PT_U2];
        // (3<   ><   >) 3'
        p += base_branch_cost * dp[st + 2][en - 1][PT_U2] * em->dangle3[stb][st1b][enb].Boltz();
        // (<   ><   >5) 5'
        p += base_branch_cost * dp[st + 1][en - 2][PT_U2] * em->dangle5[stb][en1b][enb].Boltz();
        // (.<   ><   >.) Terminal mismatch
        p += base_branch_cost * dp[st + 2][en - 2][PT_U2] *
            em->terminal[stb][st1b][en1b][enb].Boltz();

        const auto outer_coax = em->MismatchCoaxial(stb, st1b, en1b, enb);
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
              (pc.augubranch[st2b][pl1b] + em->MismatchCoaxial(pl1b, plb, st1b, st2b)).Boltz();
          // (   .(   ).) Right inner coax
          p += base_branch_cost * dp[st + 1][piv][PT_U] * dp[piv + 2][en - 2][PT_P] *
              (pc.augubranch[pr1b][en2b] + em->MismatchCoaxial(en2b, en1b, prb, pr1b)).Boltz();

          // ((   )   ) Left flush coax
          p += base_branch_cost * dp[st + 1][piv][PT_P] * dp[piv + 1][en - 1][PT_U] *
              (pc.augubranch[st1b][plb] + em->stack[stb][st1b][plb][enb]).Boltz();
          // (   (   )) Right flush coax
          p += base_branch_cost * dp[st + 1][piv][PT_U] * dp[piv + 1][en - 1][PT_P] *
              (pc.augubranch[prb][en1b] + em->stack[stb][prb][en1b][enb]).Boltz();
        }

        dp[st][en][PT_P] = p;
      }
      BoltzEnergy u{0};
      BoltzEnergy u2{0};
      BoltzEnergy rcoax{0};
      BoltzEnergy wc{0};
      BoltzEnergy gu{0};
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u += dp[st + 1][en][PT_U];
        u2 += dp[st + 1][en][PT_U2];
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

        // (   )<   > - U, U_WC?, U_GU?
        u2 += base00 * dp[piv + 1][en][PT_U];
        BoltzEnergy val = base00 + base00 * dp[piv + 1][en][PT_U];
        u += val;
        if (IsGuPair(stb, pb))
          gu += val;
        else
          wc += val;

        // (   )3<   > 3' - U
        val = base01 * em->dangle3[pl1b][pb][stb].Boltz();
        u += val;
        val *= dp[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // 5(   )<   > 5' - U
        val = base10 * em->dangle5[pb][stb][st1b].Boltz();
        u += val;
        val *= dp[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<   > Terminal mismatch - U
        val = base11 * em->terminal[pl1b][pb][stb][st1b].Boltz();
        u += val;
        val *= dp[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<(   ) > Left coax - U
        val = base11 * em->MismatchCoaxial(pl1b, pb, stb, st1b).Boltz();
        val = val * (dp[piv + 1][en][PT_U_WC] + dp[piv + 1][en][PT_U_GU]);
        u += val;
        u2 += val;

        // (   )<.(   ). > Right coax forward and backward
        val = base00 * dp[piv + 1][en][PT_U_RC];
        u += val;
        u2 += val;
        val = base11 * em->MismatchCoaxial(pl1b, pb, stb, st1b).Boltz();
        rcoax += val;
        rcoax += val * dp[piv + 1][en][PT_U];

        // (   )(<   ) > Flush coax - U
        val = base01 * em->stack[pl1b][pb][WcPair(pb)][stb].Boltz() * dp[piv][en][PT_U_WC];
        u += val;
        u2 += val;
        if (IsGu(pb)) {
          val = base01 * em->stack[pl1b][pb][GuPair(pb)][stb].Boltz() * dp[piv][en][PT_U_GU];
          u += val;
          u2 += val;
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
  PartitionExterior(r, *em, state);
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

      if (em->CanPair(r, en, st)) {
        BoltzEnergy p{0};
        const int ost_max = std::min(st + TWOLOOP_MAX_SZ + 2, N);
        for (int ost = st + 1; ost < ost_max; ++ost) {
          const int oen_min = std::max(en - TWOLOOP_MAX_SZ - 1 + (ost - st - 1), 0);
          for (int oen = en - 1; oen >= oen_min; --oen)
            p += pc.TwoLoop(oen, ost, en, st).Boltz() * dp[ost][oen][PT_P];
        }
        const BoltzEnergy base_branch_cost =
            (pc.augubranch[stb][enb] + em->multiloop_hack_a).Boltz();
        const Energy outer_coax =
            lspace && rspace ? em->MismatchCoaxial(stb, st1b, en1b, enb) : MAX_E;
        // Try being an exterior loop - coax cases handled in the loop after this.
        {
          const BoltzEnergy augu = em->AuGuPenalty(enb, stb).Boltz();
          const BoltzEnergy rext = ext[st + 1][PTEXT_R];
          const BoltzEnergy r1ext = lspace > 1 ? ext[st + 2][PTEXT_R] : BoltzEnergy{1};
          const BoltzEnergy lext = rspace ? ext[en - 1][PTEXT_L] : BoltzEnergy{1};
          const BoltzEnergy l1ext = rspace > 1 ? ext[en - 2][PTEXT_L] : BoltzEnergy{1};

          // |<   >)   (<   >| - Exterior loop
          p += augu * lext * rext;
          if (lspace) {
            // |<   >(   )3<   >| 3' - Exterior loop
            // lspace > 0
            p += augu * lext * r1ext * em->dangle3[stb][st1b][enb].Boltz();
            // |  >5)   (<   | 5' - Enclosing loop
            if (rspace > 1)
              p += base_branch_cost * dp[st + 1][en - 2][PT_U2] *
                  em->dangle5[stb][en1b][enb].Boltz();
          }
          if (rspace) {
            // |<   >5(   )<   >| 5' - Exterior loop
            // rspace > 0
            p += augu * l1ext * rext * em->dangle5[stb][en1b][enb].Boltz();
            // |   >)   (3<  | 3' - Enclosing loop
            if (lspace > 1)
              p += base_branch_cost * dp[st + 2][en - 1][PT_U2] *
                  em->dangle3[stb][st1b][enb].Boltz();
          }
          if (lspace && rspace) {
            // |<   >m(   )m<   >| Terminal mismatch - Exterior loop
            // lspace > 0 && rspace > 0
            p += augu * l1ext * r1ext * em->terminal[stb][st1b][en1b][enb].Boltz();
            // |   >)   (<   | - Enclosing loop
            p += base_branch_cost * dp[st + 1][en - 1][PT_U2];
          }
          // |  >m)   (m<  | Terminal mismatch - Enclosing loop
          if (lspace > 1 && rspace > 1)
            p += base_branch_cost * dp[st + 2][en - 2][PT_U2] *
                em->terminal[stb][st1b][en1b][enb].Boltz();

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
              p += dp[st + 2][pl][PT_P] * lext * rp1ext *
                  (em->AuGuPenalty(stb, enb) + em->AuGuPenalty(st2b, pl1b) +
                      em->MismatchCoaxial(pl1b, plb, st1b, st2b))
                      .Boltz();
              // |<   >)   ((   )<   >| Exterior loop - Left flush coax
              // lspace > 0 && not enclosed
              p += dp[st + 1][piv][PT_P] * lext * rp1ext *
                  (em->AuGuPenalty(stb, enb) + em->AuGuPenalty(st1b, plb) +
                      em->stack[stb][st1b][plb][enb])
                      .Boltz();

              if (rspace) {
                // |<   >.)   (.(   )<   >| Exterior loop - Left outer coax
                // lspace > 0 && rspace > 0 && not enclosed
                p += dp[st + 2][piv][PT_P] * l1ext * rp1ext *
                    (em->AuGuPenalty(stb, enb) + em->AuGuPenalty(st2b, plb) + outer_coax).Boltz();
              }
            }

            if (right_formable) {
              const BoltzEnergy lpext = piv > 0 ? ext[piv - 1][PTEXT_L] : BoltzEnergy{1};
              // |<   >.(   ).)   (<   >| Exterior loop - Right inner coax
              // rspace > 1 && not enclosed
              p += dp[pr][en - 2][PT_P] * lpext * rext *
                  (em->AuGuPenalty(stb, enb) + em->AuGuPenalty(prb, en2b) +
                      em->MismatchCoaxial(en2b, en1b, plb, prb))
                      .Boltz();
              // |<   >(   ))   (<   >| Exterior loop - Right flush coax
              // rspace > 0 && not enclosed
              p += dp[piv][en - 1][PT_P] * lpext * rext *
                  (em->AuGuPenalty(stb, enb) + em->AuGuPenalty(plb, en1b) +
                      em->stack[en1b][enb][stb][plb])
                      .Boltz();

              if (lspace) {
                // |<   >(   ).)   (.<   >| Exterior loop - Right outer coax
                // lspace > 0 && rspace > 1 && not enclosed
                p += dp[piv][en - 2][PT_P] * lpext * r1ext *
                    (em->AuGuPenalty(stb, enb) + em->AuGuPenalty(plb, en2b) + outer_coax).Boltz();
              }
            }
          }
        }

        // Enclosing loop cases.
        // Can start at st + 2 because we need to forman enclosing loop.
        // At worst the enclosing loop has to start at st + 1.
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
          // Left loop formable if not straddling, and is big enough or crosses over.
          const bool left_formable = straddling && (piv - st - 2 >= HAIRPIN_MIN_SZ || tpiv >= N);
          const bool right_formable =
              straddling && (en - piv - 3 >= HAIRPIN_MIN_SZ || tpiv < N - 1);
          const bool left_dot_formable = left_formable && tpiv != N;  // Can't split a dot.
          const bool right_dot_formable = right_formable && tpiv != N - 2;  // Can't split a dot.

          if (lspace > 1 && rspace > 1) {
            if (left_formable) {
              // |  >.)   (.(   )<  | Enclosing loop - Left outer coax
              // lspace > 1 && rspace > 1 && enclosed
              p += base_branch_cost * dp[st + 2][piv][PT_P] * dp[pr][en - 2][PT_U] *
                  (pc.augubranch[st2b][plb] + outer_coax).Boltz();
            }

            if (right_formable) {
              // |  >(   ).)   (.<  | Enclosing loop - Right outer coax
              // lspace > 1 && rspace > 1 && enclosed
              p += base_branch_cost * dp[st + 2][piv][PT_U] * dp[pr][en - 2][PT_P] *
                  (pc.augubranch[prb][en2b] + outer_coax).Boltz();
            }
          }

          if (lspace > 1 && rspace && left_dot_formable) {
            // |  >)   (.(   ).<  | Enclosing loop - Left inner coax
            // lspace > 1 && rspace > 0 && enclosed && no dot split
            p += base_branch_cost * dp[st + 2][pl][PT_P] * dp[pr][en - 1][PT_U] *
                (pc.augubranch[st2b][pl1b] + em->MismatchCoaxial(pl1b, plb, st1b, st2b)).Boltz();
          }

          if (lspace && rspace > 1 && right_dot_formable) {
            // |  >.(   ).)   (<  | Enclosing loop - Right inner coax
            // lspace > 0 && rspace > 1 && enclosed && no dot split
            p += base_branch_cost * dp[st + 1][piv][PT_U] * dp[pr1][en - 2][PT_P] *
                (pc.augubranch[pr1b][en2b] + em->MismatchCoaxial(en2b, en1b, prb, pr1b)).Boltz();
          }

          if (lspace && rspace) {
            if (left_formable) {
              // |  >)   ((   )<  | Enclosing loop - Left flush coax
              // lspace > 0 && rspace > 0 && enclosed
              p += base_branch_cost * dp[st + 1][piv][PT_P] * dp[pr][en - 1][PT_U] *
                  (pc.augubranch[st1b][plb] + em->stack[stb][st1b][plb][enb]).Boltz();
            }

            if (right_formable) {
              // |  >(   ))   (<  | Enclosing loop - Right flush coax
              // lspace > 0 && rspace > 0 && enclosed
              p += base_branch_cost * dp[st + 1][piv][PT_U] * dp[pr][en - 1][PT_P] *
                  (pc.augubranch[prb][en1b] + em->stack[en1b][enb][stb][prb]).Boltz();
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
      // Choose |st| to be unpaired, but only if we can maintain the constraint that we have
      // an enclosing loop formed.
      if (st + 1 < N) {
        u += dp[st + 1][en][PT_U];
        u2 += dp[st + 1][en][PT_U2];
      }

      for (int tpiv = st + 1; tpiv <= en + N; ++tpiv) {
        const int pl = FastMod(tpiv - 1, N);
        const int piv = FastMod(tpiv, N);
        const int pr = FastMod(tpiv + 1, N);
        const auto pb = r[piv];
        const auto pl1b = r[pl];
        const auto prb = r[pr];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const BoltzEnergy base00 = dp[st][piv][PT_P] * pc.augubranch[stb][pb].Boltz();
        const BoltzEnergy base01 = dp[st][pl][PT_P] * pc.augubranch[stb][pl1b].Boltz();

        // Must have an enclosing loop.
        const bool straddling = tpiv != N - 1;
        const bool dot_straddling = straddling && tpiv != N;
        BoltzEnergy val{0};

        if (straddling) {
          // |  >>   <(   )<  |
          val = base00 * dp[pr][en][PT_U];
          u2 += val;
          u += val;
          if (IsGuPair(stb, pb))
            gu += val;
          else
            wc += val;
          // U must cross the boundary to have the rest of it be nothing.
          if (tpiv >= N) {
            u += base00;
            if (IsGuPair(stb, pb))
              gu += base00;
            else
              wc += base00;
          }

          // |  )  >>   <(   )<(  | Flush coax
          // straddling
          val = base00 * em->stack[pb][prb][WcPair(prb)][stb].Boltz() * dp[pr][en][PT_U_WC];
          u += val;
          u2 += val;
          if (IsGu(prb)) {
            val = base00 * em->stack[pb][prb][GuPair(prb)][stb].Boltz() * dp[pr][en][PT_U_GU];
            u += val;
            u2 += val;
          }
          // |  ).  >>   <(   )<.(   | Right coax forward
          // straddling
          val = base00 * dp[pr][en][PT_U_RC];
          u += val;
          u2 += val;
        }

        if (dot_straddling) {
          // |  >>   <(   )3<  | 3'
          val = base01 * em->dangle3[pl1b][pb][stb].Boltz();
          if (tpiv >= N) u += val;
          val *= dp[pr][en][PT_U];
          u += val;
          u2 += val;
        }

        if (lspace) {
          const BoltzEnergy base10 = dp[st + 1][piv][PT_P] * pc.augubranch[st1b][pb].Boltz();
          const BoltzEnergy base11 = dp[st + 1][pl][PT_P] * pc.augubranch[st1b][pl1b].Boltz();

          if (straddling) {
            // |  >>   <5(   )<  | 5'
            val = base10 * em->dangle5[pb][stb][st1b].Boltz();
            if (tpiv >= N) u += val;
            val *= dp[pr][en][PT_U];
            u += val;
            u2 += val;
          }

          if (dot_straddling) {
            // |  >>   <.(   ).<  | Terminal mismatch
            // lspace > 0 && dot_straddling
            val = base11 * em->terminal[pl1b][pb][stb][st1b].Boltz();
            if (tpiv >= N) u += val;
            val *= dp[pr][en][PT_U];
            u += val;
            u2 += val;
            // |  )>>   <.(   ).<(  | Left coax
            // lspace > 0 && dot_straddling
            val = base11 * em->MismatchCoaxial(pl1b, pb, stb, st1b).Boltz();
            val = val * (dp[pr][en][PT_U_WC] + dp[pr][en][PT_U_GU]);
            u += val;
            u2 += val;
            // |  ).  >>   <(   )<.(   | Right coax backward
            val = base11 * em->MismatchCoaxial(pl1b, pb, stb, st1b).Boltz();
            if (tpiv >= N) rcoax += val;
            rcoax += val * dp[pr][en][PT_U];
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

}  // namespace mrna::md::t04
