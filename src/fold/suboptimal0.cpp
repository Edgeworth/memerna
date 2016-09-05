#include <set>
#include "fold/suboptimal0.h"

namespace memerna {
namespace fold {

using namespace constants;
using namespace energy;

Suboptimal0::Suboptimal0(const Context& ctx_, energy_t max_energy_, int max_structures_)
    : ctx(ctx_), max_energy(max_energy_), max_structures(max_structures_) {}

std::vector<computed_t> Suboptimal0::Run() {
  const auto& arr = ctx.GetDpState();
  const auto& exterior = ctx.GetExteriorState();
  const auto& em = ctx.GetEnergyModel();
  const auto& r = ctx.GetPrimary();
  const int N = int(r.size());

  // Basic idea of suboptimal traceback is look at all possible choices from a state, and expand just one of them.
  // Fully expanding one of them means there will be no duplicates in the tree.
  // Cull the ones not inside the window or when we have more than |max_structures|.
  // We don't have to check for expanding impossible states indirectly, since they will have MAX_E, be
  // above max_energy, and be instantly culled (callers use CAP_E for no energy limit).
  verify_expr(max_structures > 0, "must request at least one structure");

  q.insert({
      {{0, -1, EXT}},
      {},
      std::vector<int>(r.size(), -1),
      std::vector<Ctd>(r.size(), CTD_NA),
      exterior[0][EXT]
  });
  while (!q.empty()) {
//    printf("q size: %zu\n", q.size());
//    for (const auto& node : q) {
//      printf("  %de %zu: ", node.energy, node.not_yet_expanded.size());
//      for (const auto& nye : node.not_yet_expanded)
//        printf("(%d, %d, %d), ", nye.st, nye.en, nye.a);
//      printf("\n");
//    }
    auto node = std::move(*q.begin());
    q.erase(q.begin());
    // Finished state.
    if (node.not_yet_expanded.empty()) {
      PruneInsert(finished, node);
      continue;
    }

    // If we found a non-finished node, but |finished| is full, and the worst in |finished| is as good as
    // our current node (which is the best in |q|), then we can exit.
    if (int(finished.size()) >= max_structures && (--finished.end())->energy <= node.energy)
      break;

    const auto& to_expand = node.not_yet_expanded.back();
    node.not_yet_expanded.pop_back();
    node.history.push_back(to_expand);  // Add to history.
    int st = to_expand.st, en = to_expand.en, a = to_expand.a;

    // Initialise - we only make small modifications to it.
    curnode = node;
    // Temporary variable to hold energy calculations.
    energy_t energy = 0;

    // Exterior loop
    if (en == -1) {
      // We try replacing what we do at (st, a) with a bunch of different cases, so we use this energy as a base.
      energy_t base_energy = node.energy - exterior[st][a];
      if (a == EXT) {
        // Base case: do nothing.
        if (st == N)
          Expand(base_energy);
        else
          // Case: No pair starting here (for EXT only)
          Expand(base_energy + exterior[st + 1][EXT], {st + 1, -1, EXT});
      }
      for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        // .   .   .   (   .   .   .   )   <   >
        //           stb  st1b   en1b  enb   rem
        const auto stb = r[st], st1b = r[st + 1], enb = r[en], en1b = r[en - 1];
        const auto base00 = arr[st][en][DP_P] + em.AuGuPenalty(stb, enb);
        const auto base01 = arr[st][en - 1][DP_P] + em.AuGuPenalty(stb, en1b);
        const auto base10 = arr[st + 1][en][DP_P] + em.AuGuPenalty(st1b, enb);
        const auto base11 = arr[st + 1][en - 1][DP_P] + em.AuGuPenalty(st1b, en1b);
        curnode = node;

        // (   ).<( * ). > Right coax backward
        if (st > 0 && a == EXT_RCOAX) {
          energy = base_energy + base01 + em.MismatchCoaxial(
              en1b, enb, r[st - 1], stb) + exterior[en + 1][EXT];
          // We don't set base_ctds here, since we already set them in the forward case.
          Expand(energy, {en + 1, -1, EXT}, {st, en - 1, DP_P});
        }

        // Cases for EXT, EXT_WC, EXT_GU.
        // (   )<   >
        // If we are at EXT then this is unused.
        energy = base_energy + base00 + exterior[en + 1][EXT];
        if (a == EXT)
          Expand(energy, {en + 1, -1, EXT}, {st, en, DP_P}, {CTD_UNUSED, st});

        // (   )<   >
        // If we are at EXT_WC or EXT_GU, the CTDs for this have already have been set from a coaxial stack.
        if ((a == EXT_WC && IsWatsonCrick(stb, enb)) || (a == EXT_GU && IsGu(stb, enb)))
          Expand(energy, {en + 1, -1, EXT}, {st, en, DP_P});

        // Everything after this is only for EXT.
        if (a != EXT) continue;

        // (   )3<   > 3'
        energy = base_energy + base01 + em.dangle3[en1b][enb][stb] + exterior[en + 1][EXT];
        Expand(energy, {en + 1, -1, EXT}, {st, en - 1, DP_P}, {CTD_3_DANGLE, st});

        // 5(   )<   > 5'
        energy = base_energy + base10 + em.dangle5[enb][stb][st1b] + exterior[en + 1][EXT];
        Expand(energy, {en + 1, -1, EXT}, {st + 1, en, DP_P}, {CTD_5_DANGLE, st + 1});

        // .(   ).<   > Terminal mismatch
        energy = base_energy + base11 + em.terminal[en1b][enb][stb][st1b] + exterior[en + 1][EXT];
        Expand(energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P}, {CTD_TERMINAL_MISMATCH, st + 1});

        if (en < N - 1) {
          // .(   ).<(   ) > Left coax
          energy = base_energy + base11 + em.MismatchCoaxial(en1b, enb, stb, st1b);
          Expand(energy + exterior[en + 1][EXT_GU], {en + 1, -1, EXT_GU}, {st + 1, en - 1, DP_P},
              {CTD_LEFT_MISMATCH_COAX_WITH_PREV, en + 1}, {CTD_LEFT_MISMATCH_COAX_WITH_NEXT, st + 1});
          Expand(energy + exterior[en + 1][EXT_WC], {en + 1, -1, EXT_WC}, {st + 1, en - 1, DP_P},
              {CTD_LEFT_MISMATCH_COAX_WITH_PREV, en + 1}, {CTD_LEFT_MISMATCH_COAX_WITH_NEXT, st});

          // (   ).<(   ). > Right coax forward
          energy = base_energy + base01 + exterior[en + 1][EXT_RCOAX];
          Expand(energy, {en + 1, -1, EXT_RCOAX}, {st, en - 1, DP_P},
              {CTD_RIGHT_MISMATCH_COAX_WITH_PREV, en + 1}, {CTD_RIGHT_MISMATCH_COAX_WITH_NEXT, st});

          // (   )<(   ) > Flush coax
          const auto enrb = r[en + 1];
          energy = base_energy + base00 + em.stack[enb][enrb][enrb ^ 3][stb] + exterior[en + 1][EXT_WC];
          Expand(energy, {en + 1, -1, EXT_WC}, {st, en, DP_P},
              {CTD_FLUSH_COAX_WITH_PREV, en + 1}, {CTD_FLUSH_COAX_WITH_NEXT, st});

          if (enrb == G || enrb == U) {
            energy = base_energy + base00 + em.stack[enb][enrb][enrb ^ 1][stb] + exterior[en + 1][EXT_GU];
            Expand(energy, {en + 1, -1, EXT_GU}, {st, en, DP_P},
                {CTD_FLUSH_COAX_WITH_PREV, en + 1}, {CTD_FLUSH_COAX_WITH_NEXT, st});
          }
        }
      }
      // Finished exterior loop, don't do anymore.
      continue;
    }

    // Subtract the minimum energy of the contribution at this node.
    energy_t base_energy = node.energy - arr[st][en][a];
    // Declare the usual base aliases.
    const auto stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];

    // Normal stuff
    if (a == DP_P) {
      curnode.p[st] = en;
      curnode.p[en] = st;

      // Two loops.
      int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
      for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
        for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
          energy = base_energy + em.TwoLoop(r, st, en, ist, ien) + arr[ist][ien][DP_P];
          Expand(energy, {ist, ien, DP_P});
        }
      }

      // Hairpin loop
      energy = base_energy + em.Hairpin(r, st, en);
      Expand(energy);

      auto base_and_branch = base_energy + em.AuGuPenalty(stb, enb) + em.multiloop_hack_a + em.multiloop_hack_b;
      // (<   ><    >)
      energy = base_and_branch + arr[st + 1][en - 1][DP_U2];
      Expand(energy, {st + 1, en - 1, DP_U2}, {CTD_UNUSED, en});
      // (3<   ><   >) 3'
      energy = base_and_branch + arr[st + 2][en - 1][DP_U2] + em.dangle3[stb][st1b][enb];
      Expand(energy, {st + 2, en - 1, DP_U2}, {CTD_3_DANGLE, en});
      // (<   ><   >5) 5'
      energy = base_and_branch + arr[st + 1][en - 2][DP_U2] + em.dangle5[stb][en1b][enb];
      Expand(energy, {st + 1, en - 2, DP_U2}, {CTD_5_DANGLE, en});
      // (.<   ><   >.) Terminal mismatch
      energy = base_and_branch + arr[st + 2][en - 2][DP_U2] + em.terminal[stb][st1b][en1b][enb];
      Expand(energy, {st + 2, en - 2, DP_U2}, {CTD_TERMINAL_MISMATCH, en});

      for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
        base_t pl1b = r[piv - 1], plb = r[piv], prb = r[piv + 1], pr1b = r[piv + 2];

        // (.(   )   .) Left outer coax - P
        auto outer_coax = em.MismatchCoaxial(stb, st1b, en1b, enb);
        energy = base_and_branch + arr[st + 2][piv][DP_P] + em.multiloop_hack_b +
            em.AuGuPenalty(st2b, plb) + arr[piv + 1][en - 2][DP_U] + outer_coax;
        Expand(energy, {st + 2, piv, DP_P}, {piv + 1, en - 2, DP_U},
            {CTD_LEFT_MISMATCH_COAX_WITH_PREV, st + 1}, {CTD_LEFT_MISMATCH_COAX_WITH_NEXT, en});

        // (.   (   ).) Right outer coax
        energy = base_and_branch + arr[st + 2][piv][DP_U] + em.multiloop_hack_b +
            em.AuGuPenalty(prb, en2b) + arr[piv + 1][en - 2][DP_P] + outer_coax;
        Expand(energy, {st + 2, piv, DP_U}, {piv + 1, en - 2, DP_P},
            {CTD_RIGHT_MISMATCH_COAX_WITH_NEXT, piv + 1}, {CTD_RIGHT_MISMATCH_COAX_WITH_PREV, en});

        // (.(   ).   ) Left right coax
        energy = base_and_branch + arr[st + 2][piv - 1][DP_P] + em.multiloop_hack_b +
            em.AuGuPenalty(st2b, pl1b) + arr[piv + 1][en - 1][DP_U] +
            em.MismatchCoaxial(pl1b, plb, st1b, st2b);
        Expand(energy, {st + 2, piv - 1, DP_P}, {piv + 1, en - 1, DP_U},
            {CTD_RIGHT_MISMATCH_COAX_WITH_PREV, st + 2}, {CTD_RIGHT_MISMATCH_COAX_WITH_NEXT, en});

        // (   .(   ).) Right left coax
        energy = base_and_branch + arr[st + 1][piv][DP_U] + em.multiloop_hack_b +
            em.AuGuPenalty(pr1b, en2b) + arr[piv + 2][en - 2][DP_P] +
            em.MismatchCoaxial(en2b, en1b, prb, pr1b);
        Expand(energy, {st + 1, piv, DP_U}, {piv + 2, en - 2, DP_P},
            {CTD_LEFT_MISMATCH_COAX_WITH_NEXT, piv + 2}, {CTD_LEFT_MISMATCH_COAX_WITH_PREV, en});

        // ((   )   ) Left flush coax
        energy = base_and_branch + arr[st + 1][piv][DP_P] +
            em.multiloop_hack_b + em.AuGuPenalty(st1b, plb) +
            arr[piv + 1][en - 1][DP_U] + em.stack[stb][st1b][plb][enb];
        Expand(energy, {st + 1, piv, DP_P}, {piv + 1, en - 1, DP_U},
            {CTD_FLUSH_COAX_WITH_PREV, st + 1}, {CTD_FLUSH_COAX_WITH_NEXT, en});

        // (   (   )) Right flush coax
        energy = base_and_branch + arr[st + 1][piv][DP_U] +
            em.multiloop_hack_b + em.AuGuPenalty(prb, en1b) +
            arr[piv + 1][en - 1][DP_P] + em.stack[stb][prb][en1b][enb];
        Expand(energy, {st + 1, piv, DP_U}, {piv + 1, en - 1, DP_P},
            {CTD_FLUSH_COAX_WITH_NEXT, piv + 1}, {CTD_FLUSH_COAX_WITH_PREV, en});
      }
    } else {
      // Left unpaired. Either DP_U or DP_U2.
      if (st + 1 < en && (a == DP_U || a == DP_U2)) {
        energy = base_energy + arr[st + 1][en][a];
        Expand(energy, {st + 1, en, a});
      }

      // Pair here.
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        auto pb = r[piv], pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[st][piv][DP_P] + em.AuGuPenalty(stb, pb) + em.multiloop_hack_b;
        auto base01 = arr[st][piv - 1][DP_P] + em.AuGuPenalty(stb, pl1b) + em.multiloop_hack_b;
        auto base10 = arr[st + 1][piv][DP_P] + em.AuGuPenalty(st1b, pb) + em.multiloop_hack_b;
        auto base11 = arr[st + 1][piv - 1][DP_P] + em.AuGuPenalty(st1b, pl1b) + em.multiloop_hack_b;

        // Check a == U_RCOAX:
        // (   ).<( ** ). > Right coax backward
        if (a == DP_U_RCOAX) {
          if (st > 0) {
            energy = base_energy + base01 + em.MismatchCoaxial(pl1b, pb, r[st - 1], stb);
            // Our ctds will have already been set by now.
            Expand(energy, {st, piv - 1, DP_P});
            Expand(energy + arr[piv + 1][en][DP_U], {st, piv - 1, DP_P}, {piv + 1, en, DP_U});
          }
          continue;
        }
        // From here on, a must be U, U2, U_WC, or U_GU.

        // (   )<   > - U, U2, U_WC?, U_GU?
        energy = base_energy + base00;
        if (a == DP_U) {
          Expand(energy, {st, piv, DP_P}, {CTD_UNUSED, st});
          Expand(energy + arr[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U}, {CTD_UNUSED, st});
        }
        if (a == DP_U2)
          Expand(energy + arr[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U}, {CTD_UNUSED, st});
        if (a == DP_U_WC || a == DP_U_GU) {
          // Make sure we don't form any branches that are not the right type of pair.
          if ((a == DP_U_WC && IsWatsonCrick(stb, pb)) || (a == DP_U_GU && IsGu(stb, pb))) {
            Expand(energy, {st, piv, DP_P});
            Expand(energy + arr[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U});
          }
          continue;
        }
        // From here on, a must be U or U2.
        assert(a == DP_U || a == DP_U2);

        // (   )3<   > 3' - U, U2
        energy = base_energy + base01 + em.dangle3[pl1b][pb][stb];
        // Can only let the rest be unpaired if we only need one branch, i.e. DP_U not DP_U2.
        if (a == DP_U) Expand(energy, {st, piv - 1, DP_P}, {CTD_3_DANGLE, st});
        Expand(energy + arr[piv + 1][en][DP_U], {st, piv - 1, DP_P}, {piv + 1, en, DP_U}, {CTD_3_DANGLE, st});

        // 5(   )<   > 5' - U, U2
        energy = base_energy + base10 + em.dangle5[pb][stb][st1b];
        if (a == DP_U) Expand(energy, {st + 1, piv, DP_P}, {CTD_3_DANGLE, st + 1});
        Expand(energy + arr[piv + 1][en][DP_U], {st + 1, piv, DP_P}, {piv + 1, en, DP_U}, {CTD_3_DANGLE, st + 1});

        // .(   ).<   > Terminal mismatch - U, U2
        energy = base_energy + base11 + em.terminal[pl1b][pb][stb][st1b];
        if (a == DP_U) Expand(energy, {st + 1, piv - 1, DP_P}, {CTD_TERMINAL_MISMATCH, st + 1});
        Expand(energy + arr[piv + 1][en][DP_U], {st + 1, piv - 1, DP_P},
            {piv + 1, en, DP_U}, {CTD_TERMINAL_MISMATCH, st + 1});

        // .(   ).<(   ) > Left coax - U, U2
        energy = base_energy + base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b);
        Expand(energy + arr[piv + 1][en][DP_U_WC], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U_WC},
            {CTD_LEFT_MISMATCH_COAX_WITH_NEXT, st + 1}, {CTD_LEFT_MISMATCH_COAX_WITH_PREV, piv + 1});
        Expand(energy + arr[piv + 1][en][DP_U_GU], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U_GU},
            {CTD_LEFT_MISMATCH_COAX_WITH_NEXT, st + 1}, {CTD_LEFT_MISMATCH_COAX_WITH_PREV, piv + 1});

        // (   ).<(   ). > Right coax forward - U, U2
        energy = base_energy + base01 + arr[piv + 1][en][DP_U_RCOAX];
        Expand(energy, {st, piv - 1, DP_P}, {piv + 1, en, DP_U_RCOAX},
            {CTD_RIGHT_MISMATCH_COAX_WITH_NEXT, st}, {CTD_RIGHT_MISMATCH_COAX_WITH_PREV, piv + 1});

        // There has to be remaining bases to even have a chance at these cases.
        if (piv < en) {
          auto pr1b = r[piv + 1];
          // (   )<(   ) > Flush coax - U, U2
          energy = base_energy + base00 + em.stack[pb][pr1b][pr1b ^ 3][stb] + arr[piv + 1][en][DP_U_WC];
          Expand(energy, {st, piv, DP_P}, {piv + 1, en, DP_U_WC},
              {CTD_FLUSH_COAX_WITH_NEXT, st}, {CTD_FLUSH_COAX_WITH_PREV, piv + 1});

          if (pr1b == G || pr1b == U) {
            energy = base_energy + base00 + em.stack[pb][pr1b][pr1b ^ 1][stb] + arr[piv + 1][en][DP_U_GU];
            Expand(energy, {st, piv, DP_P}, {piv + 1, en, DP_U_GU},
                {CTD_FLUSH_COAX_WITH_NEXT, st}, {CTD_FLUSH_COAX_WITH_PREV, piv + 1});
          }
        }
      }
    }
  }

  std::vector<computed_t> ret;
  for (const auto& struc : finished) {
    assert(struc.not_yet_expanded.empty());
    ret.push_back({{r, struc.p}, struc.base_ctds, struc.energy});
  }
  return ret;
}

}
}
