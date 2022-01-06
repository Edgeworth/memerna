// Copyright 2016 Eliot Courtney.
#include "compute/energy/globals.h"
#include "compute/partition/globals.h"
#include "compute/partition/partition.h"

namespace mrna::partition::internal {

using energy::Boltzmann;
using energy::gem;
using energy::gpc;

void Partition0() {
  const int N = static_cast<int>(gr.size());
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const Base stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en],
                 en1b = gr[en - 1], en2b = gr[en - 2];

      // if (CanPair(stb, enb)) {  // TODO lonely pairs?
      if (energy::ViableFoldingPair(st, en)) {
        PEnergy p{0};
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist)  // TODO lyngso's ?
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien)
            p += Boltzmann(energy::FastTwoLoop(st, en, ist, ien)) * gpt[ist][ien][PT_P];
        // Hairpin loops.
        p += Boltzmann(gem.Hairpin(gr, st, en));
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const PEnergy base_branch_cost = Boltzmann(gpc.augubranch[stb][enb] + gem.multiloop_hack_a);

        // (<   ><   >)
        p += base_branch_cost * gpt[st + 1][en - 1][PT_U2];
        // (3<   ><   >) 3'
        p += base_branch_cost * gpt[st + 2][en - 1][PT_U2] * Boltzmann(gem.dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        p += base_branch_cost * gpt[st + 1][en - 2][PT_U2] * Boltzmann(gem.dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        p += base_branch_cost * gpt[st + 2][en - 2][PT_U2] *
            Boltzmann(gem.terminal[stb][st1b][en1b][enb]);

        const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);
        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          Base pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1], pr1b = gr[piv + 2];

          // (.(   )   .) Left outer coax - P
          p += base_branch_cost * gpt[st + 2][piv][PT_P] * gpt[piv + 1][en - 2][PT_U] *
              Boltzmann(gpc.augubranch[st2b][plb] + outer_coax);
          // (.   (   ).) Right outer coax
          p += base_branch_cost * gpt[st + 2][piv][PT_U] * gpt[piv + 1][en - 2][PT_P] *
              Boltzmann(gpc.augubranch[prb][en2b] + outer_coax);

          // (.(   ).   ) Left right coax
          p += base_branch_cost * gpt[st + 2][piv - 1][PT_P] * gpt[piv + 1][en - 1][PT_U] *
              Boltzmann(gpc.augubranch[st2b][pl1b] + gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[piv + 2][en - 2][PT_P] *
              Boltzmann(gpc.augubranch[pr1b][en2b] + gem.MismatchCoaxial(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          p += base_branch_cost * gpt[st + 1][piv][PT_P] * gpt[piv + 1][en - 1][PT_U] *
              Boltzmann(gpc.augubranch[st1b][plb] + gem.stack[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[piv + 1][en - 1][PT_P] *
              Boltzmann(gpc.augubranch[prb][en1b] + gem.stack[stb][prb][en1b][enb]);
        }

        gpt[st][en][PT_P] = p;
      }
      PEnergy u{0}, u2{0}, rcoax{0}, wc{0}, gu{0};
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u += gpt[st + 1][en][PT_U];
        u2 += gpt[st + 1][en][PT_U2];
      }

      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = gr[piv], pl1b = gr[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const PEnergy base00 = gpt[st][piv][PT_P] * Boltzmann(gpc.augubranch[stb][pb]);
        const PEnergy base01 = gpt[st][piv - 1][PT_P] * Boltzmann(gpc.augubranch[stb][pl1b]);
        const PEnergy base10 = gpt[st + 1][piv][PT_P] * Boltzmann(gpc.augubranch[st1b][pb]);
        const PEnergy base11 = gpt[st + 1][piv - 1][PT_P] * Boltzmann(gpc.augubranch[st1b][pl1b]);

        // (   )<   > - U, U_WC?, U_GU?
        u2 += base00 * gpt[piv + 1][en][PT_U];
        PEnergy val = base00 + base00 * gpt[piv + 1][en][PT_U];
        u += val;
        if (IsGu(stb, pb))
          gu += val;
        else
          wc += val;

        // (   )3<   > 3' - U
        val = base01 * Boltzmann(gem.dangle3[pl1b][pb][stb]);
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // 5(   )<   > 5' - U
        val = base10 * Boltzmann(gem.dangle5[pb][stb][st1b]);
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<   > Terminal mismatch - U
        val = base11 * Boltzmann(gem.terminal[pl1b][pb][stb][st1b]);
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<(   ) > Left coax - U
        val = base11 * Boltzmann(gem.MismatchCoaxial(pl1b, pb, stb, st1b));
        val = val * (gpt[piv + 1][en][PT_U_WC] + gpt[piv + 1][en][PT_U_GU]);
        u += val;
        u2 += val;

        // (   )<.(   ). > Right coax forward and backward
        val = base00 * gpt[piv + 1][en][PT_U_RCOAX];
        u += val;
        u2 += val;
        val = base11 * Boltzmann(gem.MismatchCoaxial(pl1b, pb, stb, st1b));
        rcoax += val;
        rcoax += val * gpt[piv + 1][en][PT_U];

        // (   )(<   ) > Flush coax - U
        val = base01 * Boltzmann(gem.stack[pl1b][pb][pb ^ 3][stb]) * gpt[piv][en][PT_U_WC];
        u += val;
        u2 += val;
        if (pb == G || pb == U) {
          val = base01 * Boltzmann(gem.stack[pl1b][pb][pb ^ 1][stb]) * gpt[piv][en][PT_U_GU];
          u += val;
          u2 += val;
        }
      }
      gpt[st][en][PT_U] = u;
      gpt[st][en][PT_U2] = u2;
      gpt[st][en][PT_U_WC] = wc;
      gpt[st][en][PT_U_GU] = gu;
      gpt[st][en][PT_U_RCOAX] = rcoax;
    }
  }

  // Compute the exterior tables.
  Exterior();

  // Fill the left triangle.
  // The meaning of the tables changes here:
  // U, U2, WC, GU, RCOAX: any table index with en < st must have a loop enclosing (en, st)
  for (int st = N - 1; st >= 0; --st) {
    for (int en = 0; en < st; ++en) {
      //        ..)...(..
      // rspace  en   st  lspace
      const int lspace = N - st - 1, rspace = en;
      const Base stb = gr[st], st1b = lspace ? gr[st + 1] : Base(-1),
                 st2b = lspace > 1 ? gr[st + 2] : Base(-1), enb = gr[en],
                 en1b = rspace ? gr[en - 1] : Base(-1), en2b = rspace > 1 ? gr[en - 2] : Base(-1);

      // if (CanPair(enb, stb)) {  // TODO lonely pairs?
      if (energy::ViableFoldingPair(en, st)) {
        PEnergy p{0};
        const int ost_max = std::min(st + TWOLOOP_MAX_SZ + 2, N);
        for (int ost = st + 1; ost < ost_max; ++ost) {
          const int oen_min = std::max(en - TWOLOOP_MAX_SZ - 1 + (ost - st - 1), 0);
          for (int oen = en - 1; oen >= oen_min; --oen)
            p += Boltzmann(energy::FastTwoLoop(oen, ost, en, st)) * gpt[ost][oen][PT_P];
        }
        const PEnergy base_branch_cost = Boltzmann(gpc.augubranch[stb][enb] + gem.multiloop_hack_a);
        const Energy outer_coax =
            lspace && rspace ? gem.MismatchCoaxial(stb, st1b, en1b, enb) : MAX_E;
        // Try being an exterior loop - coax cases handled in the loop after this.
        {
          const PEnergy augu = Boltzmann(gem.AuGuPenalty(enb, stb));
          const PEnergy rext = gptext[st + 1][PTEXT_R];
          const PEnergy r1ext = lspace > 1 ? gptext[st + 2][PTEXT_R] : PEnergy{1};
          const PEnergy lext = rspace ? gptext[en - 1][PTEXT_L] : PEnergy{1};
          const PEnergy l1ext = rspace > 1 ? gptext[en - 2][PTEXT_L] : PEnergy{1};

          // |<   >)   (<   >| - Exterior loop
          p += augu * lext * rext;
          if (lspace) {
            // |<   >(   )3<   >| 3' - Exterior loop
            // lspace > 0
            p += augu * lext * r1ext * Boltzmann(gem.dangle3[stb][st1b][enb]);
            // |  >5)   (<   | 5' - Enclosing loop
            if (rspace > 1)
              p += base_branch_cost * gpt[st + 1][en - 2][PT_U2] *
                  Boltzmann(gem.dangle5[stb][en1b][enb]);
          }
          if (rspace) {
            // |<   >5(   )<   >| 5' - Exterior loop
            // rspace > 0
            p += augu * l1ext * rext * Boltzmann(gem.dangle5[stb][en1b][enb]);
            // |   >)   (3<  | 3' - Enclosing loop
            if (lspace > 1)
              p += base_branch_cost * gpt[st + 2][en - 1][PT_U2] *
                  Boltzmann(gem.dangle3[stb][st1b][enb]);
          }
          if (lspace && rspace) {
            // |<   >m(   )m<   >| Terminal mismatch - Exterior loop
            // lspace > 0 && rspace > 0
            p += augu * l1ext * r1ext * Boltzmann(gem.terminal[stb][st1b][en1b][enb]);
            // |   >)   (<   | - Enclosing loop
            p += base_branch_cost * gpt[st + 1][en - 1][PT_U2];
          }
          // |  >m)   (m<  | Terminal mismatch - Enclosing loop
          if (lspace > 1 && rspace > 1)
            p += base_branch_cost * gpt[st + 2][en - 2][PT_U2] *
                Boltzmann(gem.terminal[stb][st1b][en1b][enb]);

          const int limit = en + N;
          for (int tpiv = st; tpiv <= limit; ++tpiv) {
            const int pl = FastMod(tpiv - 1, N), piv = FastMod(tpiv, N), pr = FastMod(tpiv + 1, N);
            Base pl1b = gr[pl], plb = gr[piv], prb = gr[pr];
            const bool left_formable = tpiv < N && tpiv - st - 2 >= HAIRPIN_MIN_SZ;
            const bool right_formable = tpiv >= N && en - piv - 2 >= HAIRPIN_MIN_SZ;

            if (left_formable) {
              const PEnergy rp1ext = gptext[piv + 1][PTEXT_R];
              // |<   >)   (.(   ).<   >| Exterior loop - Left right coax
              // lspace > 1 && not enclosed
              p += gpt[st + 2][pl][PT_P] * lext * rp1ext *
                  Boltzmann(gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(st2b, pl1b) +
                      gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
              // |<   >)   ((   )<   >| Exterior loop - Left flush coax
              // lspace > 0 && not enclosed
              p += gpt[st + 1][piv][PT_P] * lext * rp1ext *
                  Boltzmann(gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(st1b, plb) +
                      gem.stack[stb][st1b][plb][enb]);

              if (rspace) {
                // |<   >.)   (.(   )<   >| Exterior loop - Left outer coax
                // lspace > 0 && rspace > 0 && not enclosed
                p += gpt[st + 2][piv][PT_P] * l1ext * rp1ext *
                    Boltzmann(gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(st2b, plb) + outer_coax);
              }
            }

            if (right_formable) {
              const PEnergy lpext = piv > 0 ? gptext[piv - 1][PTEXT_L] : PEnergy{1};
              // |<   >.(   ).)   (<   >| Exterior loop - Right left coax
              // rspace > 1 && not enclosed
              p += gpt[pr][en - 2][PT_P] * lpext * rext *
                  Boltzmann(gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(prb, en2b) +
                      gem.MismatchCoaxial(en2b, en1b, plb, prb));
              // |<   >(   ))   (<   >| Exterior loop - Right flush coax
              // rspace > 0 && not enclosed
              p += gpt[piv][en - 1][PT_P] * lpext * rext *
                  Boltzmann(gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(plb, en1b) +
                      gem.stack[en1b][enb][stb][plb]);

              if (lspace) {
                // |<   >(   ).)   (.<   >| Exterior loop - Right outer coax
                // lspace > 0 && rspace > 1 && not enclosed
                p += gpt[piv][en - 2][PT_P] * lpext * r1ext *
                    Boltzmann(gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(plb, en2b) + outer_coax);
              }
            }
          }
        }

        // Enclosing loop cases.
        // Can start at st + 2 because we need to forman enclosing loop.
        // At worst the enclosing loop has to start at st + 1.
        const int limit = en + N;
        for (int tpiv = st + 2; tpiv < limit; ++tpiv) {
          const int pl = FastMod(tpiv - 1, N), piv = FastMod(tpiv, N), pr = FastMod(tpiv + 1, N),
                    pr1 = FastMod(tpiv + 2, N);
          // Left block is: [st, piv], Right block is: [piv + 1, en].
          Base pl1b = gr[pl], plb = gr[piv], prb = gr[pr], pr1b = gr[pr1];
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
              p += base_branch_cost * gpt[st + 2][piv][PT_P] * gpt[pr][en - 2][PT_U] *
                  Boltzmann(gpc.augubranch[st2b][plb] + outer_coax);
            }

            if (right_formable) {
              // |  >(   ).)   (.<  | Enclosing loop - Right outer coax
              // lspace > 1 && rspace > 1 && enclosed
              p += base_branch_cost * gpt[st + 2][piv][PT_U] * gpt[pr][en - 2][PT_P] *
                  Boltzmann(gpc.augubranch[prb][en2b] + outer_coax);
            }
          }

          if (lspace > 1 && rspace && left_dot_formable) {
            // |  >)   (.(   ).<  | Enclosing loop - Left right coax
            // lspace > 1 && rspace > 0 && enclosed && no dot split
            p += base_branch_cost * gpt[st + 2][pl][PT_P] * gpt[pr][en - 1][PT_U] *
                Boltzmann(gpc.augubranch[st2b][pl1b] + gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
          }

          if (lspace && rspace > 1 && right_dot_formable) {
            // |  >.(   ).)   (<  | Enclosing loop - Right left coax
            // lspace > 0 && rspace > 1 && enclosed && no dot split
            p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[pr1][en - 2][PT_P] *
                Boltzmann(gpc.augubranch[pr1b][en2b] + gem.MismatchCoaxial(en2b, en1b, prb, pr1b));
          }

          if (lspace && rspace) {
            if (left_formable) {
              // |  >)   ((   )<  | Enclosing loop - Left flush coax
              // lspace > 0 && rspace > 0 && enclosed
              p += base_branch_cost * gpt[st + 1][piv][PT_P] * gpt[pr][en - 1][PT_U] *
                  Boltzmann(gpc.augubranch[st1b][plb] + gem.stack[stb][st1b][plb][enb]);
            }

            if (right_formable) {
              // |  >(   ))   (<  | Enclosing loop - Right flush coax
              // lspace > 0 && rspace > 0 && enclosed
              p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[pr][en - 1][PT_P] *
                  Boltzmann(gpc.augubranch[prb][en1b] + gem.stack[en1b][enb][stb][prb]);
            }
          }
        }
        gpt[st][en][PT_P] = p;
      }
      PEnergy u{0}, u2{0}, rcoax{0}, wc{0}, gu{0};
      // Update unpaired.
      // Choose |st| to be unpaired, but only if we can maintain the constraint that we have
      // an enclosing loop formed.
      if (st + 1 < N) {
        u += gpt[st + 1][en][PT_U];
        u2 += gpt[st + 1][en][PT_U2];
      }

      for (int tpiv = st + 1; tpiv <= en + N; ++tpiv) {
        const int pl = FastMod(tpiv - 1, N), piv = FastMod(tpiv, N), pr = FastMod(tpiv + 1, N);
        const auto pb = gr[piv], pl1b = gr[pl], prb = gr[pr];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const PEnergy base00 = gpt[st][piv][PT_P] * Boltzmann(gpc.augubranch[stb][pb]);
        const PEnergy base01 = gpt[st][pl][PT_P] * Boltzmann(gpc.augubranch[stb][pl1b]);

        // Must have an enclosing loop.
        const bool straddling = tpiv != N - 1;
        const bool dot_straddling = straddling && tpiv != N;
        PEnergy val{0};

        if (straddling) {
          // |  >>   <(   )<  |
          val = base00 * gpt[pr][en][PT_U];
          u2 += val;
          u += val;
          if (IsGu(stb, pb))
            gu += val;
          else
            wc += val;
          // U must cross the boundary to have the rest of it be nothing.
          if (tpiv >= N) {
            u += base00;
            if (IsGu(stb, pb))
              gu += base00;
            else
              wc += base00;
          }

          // |  )  >>   <(   )<(  | Flush coax
          // straddling
          val = base00 * Boltzmann(gem.stack[pb][prb][prb ^ 3][stb]) * gpt[pr][en][PT_U_WC];
          u += val;
          u2 += val;
          if (prb == G || prb == U) {
            val = base00 * Boltzmann(gem.stack[pb][prb][prb ^ 1][stb]) * gpt[pr][en][PT_U_GU];
            u += val;
            u2 += val;
          }
          // |  ).  >>   <(   )<.(   | Right coax forward
          // straddling
          val = base00 * gpt[pr][en][PT_U_RCOAX];
          u += val;
          u2 += val;
        }

        if (dot_straddling) {
          // |  >>   <(   )3<  | 3'
          val = base01 * Boltzmann(gem.dangle3[pl1b][pb][stb]);
          if (tpiv >= N) u += val;
          val *= gpt[pr][en][PT_U];
          u += val;
          u2 += val;
        }

        if (lspace) {
          const PEnergy base10 = gpt[st + 1][piv][PT_P] * Boltzmann(gpc.augubranch[st1b][pb]);
          const PEnergy base11 = gpt[st + 1][pl][PT_P] * Boltzmann(gpc.augubranch[st1b][pl1b]);

          if (straddling) {
            // |  >>   <5(   )<  | 5'
            val = base10 * Boltzmann(gem.dangle5[pb][stb][st1b]);
            if (tpiv >= N) u += val;
            val *= gpt[pr][en][PT_U];
            u += val;
            u2 += val;
          }

          if (dot_straddling) {
            // |  >>   <.(   ).<  | Terminal mismatch
            // lspace > 0 && dot_straddling
            val = base11 * Boltzmann(gem.terminal[pl1b][pb][stb][st1b]);
            if (tpiv >= N) u += val;
            val *= gpt[pr][en][PT_U];
            u += val;
            u2 += val;
            // |  )>>   <.(   ).<(  | Left coax
            // lspace > 0 && dot_straddling
            val = base11 * Boltzmann(gem.MismatchCoaxial(pl1b, pb, stb, st1b));
            val = val * (gpt[pr][en][PT_U_WC] + gpt[pr][en][PT_U_GU]);
            u += val;
            u2 += val;
            // |  ).  >>   <(   )<.(   | Right coax backward
            val = base11 * Boltzmann(gem.MismatchCoaxial(pl1b, pb, stb, st1b));
            if (tpiv >= N) rcoax += val;
            rcoax += val * gpt[pr][en][PT_U];
          }
        }
      }

      gpt[st][en][PT_U] = u;
      gpt[st][en][PT_U2] = u2;
      gpt[st][en][PT_U_WC] = wc;
      gpt[st][en][PT_U_GU] = gu;
      gpt[st][en][PT_U_RCOAX] = rcoax;
    }
  }
}

}  // namespace mrna::partition::internal
