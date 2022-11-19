// Copyright 2016 E.
#include "compute/subopt/subopt0.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "compute/subopt/subopt.h"
#include "compute/traceback/traceback.h"
#include "model/base.h"
#include "model/secondary.h"
#include "util/array.h"
#include "util/error.h"

namespace mrna::subopt {

Suboptimal0::Suboptimal0(
    Primary r, energy::EnergyModelPtr em, DpArray dp, ExtArray ext, SuboptCfg cfg)
    : r_(std::move(r)), em_(std::move(em)), dp_(std::move(dp)), ext_(std::move(ext)), cfg_(cfg) {}

int Suboptimal0::Run(const SuboptCallback& fn) {
  const int N = static_cast<int>(r_.size());
  verify(N < std::numeric_limits<int16_t>::max(), "RNA too long for suboptimal folding");

  // Basic idea of suboptimal traceback is look at all possible choices from a state, and expand
  // just one of them. Fully expanding one of them means there will be no duplicates in the tree.
  // Cull the ones not inside the window or when we have more than |max_structures|.
  // We don't have to check for expanding impossible states indirectly, since they will have MAX_E,
  // be above cfg_.delta, and be instantly culled (callers use CAP_E for no energy limit).
  q_.insert({.not_yet_expanded = {{0, -1, EXT}},
      .history{},
      .res = SuboptResult(ext_[0][EXT], tb::TracebackResult(Secondary(N), Ctds(N)))});
  Node node;
  while (!q_.empty()) {
    node = std::move(q_.extract(q_.begin()).value());
    // Finished state.
    if (node.not_yet_expanded.empty()) {
      PruneInsert(node, &finished_);
      continue;
    }

    // If we found a non-finished node, but |finished| is full, and the worst in |finished| is
    // as good as our current node (which is the best in |q|), then we can exit.
    if (static_cast<int>(finished_.size()) >= cfg_.strucs &&
        (--finished_.end())->res.energy <= node.res.energy)
      break;

    auto to_expand = node.not_yet_expanded.back();
    node.not_yet_expanded.pop_back();
    node.history.push_back(to_expand);  // Add to history.
    int st = to_expand.st;
    int en = to_expand.en;
    int a = to_expand.a;

    // Initialise - we only make small modifications to it.
    curnode_ = node.copy();
    // Temporary variable to hold energy calculations.
    Energy energy = ZERO_E;

    // Exterior loop
    if (en == -1) {
      // We try replacing what we do at (st, a) with a bunch of different cases, so we use this
      // energy as a base.
      Energy base_energy = node.res.energy - ext_[st][a];
      if (a == EXT) {
        // Base case: do nothing.
        if (st == N)
          Expand(base_energy);
        else
          // Case: No pair starting here (for EXT only)
          Expand(base_energy + ext_[st + 1][EXT], {st + 1, -1, EXT});
      }
      for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        // .   .   .   (   .   .   .   )   <   >
        //           stb  st1b   en1b  enb   rem
        const auto stb = r_[st];
        const auto st1b = r_[st + 1];
        const auto enb = r_[en];
        const auto en1b = r_[en - 1];
        const auto base00 = dp_[st][en][DP_P] + em_->AuGuPenalty(stb, enb);
        const auto base01 = dp_[st][en - 1][DP_P] + em_->AuGuPenalty(stb, en1b);
        const auto base10 = dp_[st + 1][en][DP_P] + em_->AuGuPenalty(st1b, enb);
        const auto base11 = dp_[st + 1][en - 1][DP_P] + em_->AuGuPenalty(st1b, en1b);
        curnode_ = node.copy();

        // (   )<.( * ). > Right coax backward
        if (a == EXT_RC) {
          energy =
              base_energy + base11 + em_->MismatchCoaxial(en1b, enb, stb, st1b) + ext_[en + 1][EXT];
          // We don't set ctds here, since we already set them in the forward case.
          Expand(energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P});
        }

        // Cases for EXT, EXT_WC, EXT_GU.
        // (   )<   >
        // If we are at EXT then this is unused.
        energy = base_energy + base00 + ext_[en + 1][EXT];
        if (a == EXT) Expand(energy, {en + 1, -1, EXT}, {st, en, DP_P}, {st, CTD_UNUSED});

        // (   )<   >
        // If we are at EXT_WC or EXT_GU, the CTDs for this have already have been set from a
        // coaxial stack.
        if ((a == EXT_WC && IsWcPair(stb, enb)) || (a == EXT_GU && IsGuPair(stb, enb)))
          Expand(energy, {en + 1, -1, EXT}, {st, en, DP_P});

        // Everything after this is only for EXT.
        if (a != EXT) continue;

        // (   )3<   > 3'
        energy = base_energy + base01 + em_->dangle3[en1b][enb][stb] + ext_[en + 1][EXT];
        Expand(energy, {en + 1, -1, EXT}, {st, en - 1, DP_P}, {st, CTD_3_DANGLE});

        // 5(   )<   > 5'
        energy = base_energy + base10 + em_->dangle5[enb][stb][st1b] + ext_[en + 1][EXT];
        Expand(energy, {en + 1, -1, EXT}, {st + 1, en, DP_P}, {st + 1, CTD_5_DANGLE});

        // .(   ).<   > Terminal mismatch
        energy = base_energy + base11 + em_->terminal[en1b][enb][stb][st1b] + ext_[en + 1][EXT];
        Expand(energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P}, {st + 1, CTD_MISMATCH});

        // (   )(<   ) > Flush coax
        energy = base_energy + base01 + em_->stack[en1b][enb][WcPair(enb)][stb] + ext_[en][EXT_WC];
        Expand(energy, {en, -1, EXT_WC}, {st, en - 1, DP_P}, {en, CTD_FCOAX_WITH_PREV},
            {st, CTD_FCOAX_WITH_NEXT});
        if (IsGu(enb)) {
          energy =
              base_energy + base01 + em_->stack[en1b][enb][GuPair(enb)][stb] + ext_[en][EXT_GU];
          Expand(energy, {en, -1, EXT_GU}, {st, en - 1, DP_P}, {en, CTD_FCOAX_WITH_PREV},
              {st, CTD_FCOAX_WITH_NEXT});
        }

        if (en < N - 1) {
          // .(   ).<(   ) > Left coax
          energy = base_energy + base11 + em_->MismatchCoaxial(en1b, enb, stb, st1b);
          Expand(energy + ext_[en + 1][EXT_GU], {en + 1, -1, EXT_GU}, {st + 1, en - 1, DP_P},
              {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT});
          Expand(energy + ext_[en + 1][EXT_WC], {en + 1, -1, EXT_WC}, {st + 1, en - 1, DP_P},
              {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT});
        }
        if (en < N - 2) {
          // (   )<.(   ). > Right coax forward
          energy = base_energy + base00 + ext_[en + 1][EXT_RC];
          Expand(energy, {en + 1, -1, EXT_RC}, {st, en, DP_P}, {en + 2, CTD_RC_WITH_PREV},
              {st, CTD_RC_WITH_NEXT});
        }
      }
      // Finished exterior loop, don't do anymore.
      continue;
    }

    // Subtract the minimum energy of the contribution at this node.
    Energy base_energy = node.res.energy - dp_[st][en][a];
    // Declare the usual base aliases.
    const auto stb = r_[st];
    const auto st1b = r_[st + 1];
    const auto st2b = r_[st + 2];
    const auto enb = r_[en];
    const auto en1b = r_[en - 1];
    const auto en2b = r_[en - 2];

    // Normal stuff
    if (a == DP_P) {
      curnode_.res.tb.s[st] = en;
      curnode_.res.tb.s[en] = st;

      // Two loops.
      int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
      for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
        for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
          energy = base_energy + em_->TwoLoop(r_, st, en, ist, ien) + dp_[ist][ien][DP_P];
          Expand(energy, {ist, ien, DP_P});
        }
      }

      // Hairpin loop
      energy = base_energy + em_->Hairpin(r_, st, en);
      Expand(energy);

      auto base_and_branch =
          base_energy + em_->AuGuPenalty(stb, enb) + em_->multiloop_hack_a + em_->multiloop_hack_b;
      // (<   ><    >)
      energy = base_and_branch + dp_[st + 1][en - 1][DP_U2];
      Expand(energy, {st + 1, en - 1, DP_U2}, {en, CTD_UNUSED});
      // (3<   ><   >) 3'
      energy = base_and_branch + dp_[st + 2][en - 1][DP_U2] + em_->dangle3[stb][st1b][enb];
      Expand(energy, {st + 2, en - 1, DP_U2}, {en, CTD_3_DANGLE});
      // (<   ><   >5) 5'
      energy = base_and_branch + dp_[st + 1][en - 2][DP_U2] + em_->dangle5[stb][en1b][enb];
      Expand(energy, {st + 1, en - 2, DP_U2}, {en, CTD_5_DANGLE});
      // (.<   ><   >.) Terminal mismatch
      energy = base_and_branch + dp_[st + 2][en - 2][DP_U2] + em_->terminal[stb][st1b][en1b][enb];
      Expand(energy, {st + 2, en - 2, DP_U2}, {en, CTD_MISMATCH});

      for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
        Base pl1b = r_[piv - 1];
        Base plb = r_[piv];
        Base prb = r_[piv + 1];
        Base pr1b = r_[piv + 2];

        // (.(   )   .) Left outer coax - P
        auto outer_coax = em_->MismatchCoaxial(stb, st1b, en1b, enb);
        energy = base_and_branch + dp_[st + 2][piv][DP_P] + em_->multiloop_hack_b +
            em_->AuGuPenalty(st2b, plb) + dp_[piv + 1][en - 2][DP_U] + outer_coax;
        Expand(energy, {st + 2, piv, DP_P}, {piv + 1, en - 2, DP_U}, {st + 2, CTD_LCOAX_WITH_PREV},
            {en, CTD_LCOAX_WITH_NEXT});

        // (.   (   ).) Right outer coax
        energy = base_and_branch + dp_[st + 2][piv][DP_U] + em_->multiloop_hack_b +
            em_->AuGuPenalty(prb, en2b) + dp_[piv + 1][en - 2][DP_P] + outer_coax;
        Expand(energy, {st + 2, piv, DP_U}, {piv + 1, en - 2, DP_P}, {piv + 1, CTD_RC_WITH_NEXT},
            {en, CTD_RC_WITH_PREV});

        // (.(   ).   ) Left inner coax
        energy = base_and_branch + dp_[st + 2][piv - 1][DP_P] + em_->multiloop_hack_b +
            em_->AuGuPenalty(st2b, pl1b) + dp_[piv + 1][en - 1][DP_U] +
            em_->MismatchCoaxial(pl1b, plb, st1b, st2b);
        Expand(energy, {st + 2, piv - 1, DP_P}, {piv + 1, en - 1, DP_U}, {st + 2, CTD_RC_WITH_PREV},
            {en, CTD_RC_WITH_NEXT});

        // (   .(   ).) Right inner coax
        energy = base_and_branch + dp_[st + 1][piv][DP_U] + em_->multiloop_hack_b +
            em_->AuGuPenalty(pr1b, en2b) + dp_[piv + 2][en - 2][DP_P] +
            em_->MismatchCoaxial(en2b, en1b, prb, pr1b);
        Expand(energy, {st + 1, piv, DP_U}, {piv + 2, en - 2, DP_P}, {piv + 2, CTD_LCOAX_WITH_NEXT},
            {en, CTD_LCOAX_WITH_PREV});

        // ((   )   ) Left flush coax
        energy = base_and_branch + dp_[st + 1][piv][DP_P] + em_->multiloop_hack_b +
            em_->AuGuPenalty(st1b, plb) + dp_[piv + 1][en - 1][DP_U] +
            em_->stack[stb][st1b][plb][enb];
        Expand(energy, {st + 1, piv, DP_P}, {piv + 1, en - 1, DP_U}, {st + 1, CTD_FCOAX_WITH_PREV},
            {en, CTD_FCOAX_WITH_NEXT});

        // (   (   )) Right flush coax
        energy = base_and_branch + dp_[st + 1][piv][DP_U] + em_->multiloop_hack_b +
            em_->AuGuPenalty(prb, en1b) + dp_[piv + 1][en - 1][DP_P] +
            em_->stack[stb][prb][en1b][enb];
        Expand(energy, {st + 1, piv, DP_U}, {piv + 1, en - 1, DP_P}, {piv + 1, CTD_FCOAX_WITH_NEXT},
            {en, CTD_FCOAX_WITH_PREV});
      }
    } else {
      // Left unpaired. Either DP_U or DP_U2.
      if (st + 1 < en && (a == DP_U || a == DP_U2)) {
        energy = base_energy + dp_[st + 1][en][a];
        Expand(energy, {st + 1, en, a});
      }

      // Pair here.
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        auto pb = r_[piv];
        auto pl1b = r_[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the
        // right.
        auto base00 = dp_[st][piv][DP_P] + em_->AuGuPenalty(stb, pb) + em_->multiloop_hack_b;
        auto base01 = dp_[st][piv - 1][DP_P] + em_->AuGuPenalty(stb, pl1b) + em_->multiloop_hack_b;
        auto base10 = dp_[st + 1][piv][DP_P] + em_->AuGuPenalty(st1b, pb) + em_->multiloop_hack_b;
        auto base11 =
            dp_[st + 1][piv - 1][DP_P] + em_->AuGuPenalty(st1b, pl1b) + em_->multiloop_hack_b;

        // Check a == U_RC:
        // (   )<.( ** ). > Right coax backward
        if (a == DP_U_RC) {
          energy = base_energy + base11 + em_->MismatchCoaxial(pl1b, pb, stb, st1b);
          // Our ctds will have already been set by now.
          Expand(energy, {st + 1, piv - 1, DP_P});
          Expand(energy + dp_[piv + 1][en][DP_U], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U});
          continue;
        }
        // From here on, a must be U, U2, U_WC, or U_GU.

        // (   )<   > - U, U2, U_WC?, U_GU?
        energy = base_energy + base00;
        if (a == DP_U) {
          Expand(energy, {st, piv, DP_P}, {st, CTD_UNUSED});
          Expand(energy + dp_[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U},
              {st, CTD_UNUSED});
        }
        if (a == DP_U2)
          Expand(energy + dp_[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U},
              {st, CTD_UNUSED});
        if (a == DP_U_WC || a == DP_U_GU) {
          // Make sure we don't form any branches that are not the right type of pair.
          if ((a == DP_U_WC && IsWcPair(stb, pb)) || (a == DP_U_GU && IsGuPair(stb, pb))) {
            Expand(energy, {st, piv, DP_P});
            Expand(energy + dp_[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U});
          }
          continue;
        }
        // From here on, a must be U or U2.
        assert(a == DP_U || a == DP_U2);

        // (   )3<   > 3' - U, U2
        energy = base_energy + base01 + em_->dangle3[pl1b][pb][stb];
        // Can only let the rest be unpaired if we only need one branch, i.e. DP_U not DP_U2.
        if (a == DP_U) Expand(energy, {st, piv - 1, DP_P}, {st, CTD_3_DANGLE});
        Expand(energy + dp_[piv + 1][en][DP_U], {st, piv - 1, DP_P}, {piv + 1, en, DP_U},
            {st, CTD_3_DANGLE});

        // 5(   )<   > 5' - U, U2
        energy = base_energy + base10 + em_->dangle5[pb][stb][st1b];
        if (a == DP_U) Expand(energy, {st + 1, piv, DP_P}, {st + 1, CTD_5_DANGLE});
        Expand(energy + dp_[piv + 1][en][DP_U], {st + 1, piv, DP_P}, {piv + 1, en, DP_U},
            {st + 1, CTD_5_DANGLE});

        // .(   ).<   > Terminal mismatch - U, U2
        energy = base_energy + base11 + em_->terminal[pl1b][pb][stb][st1b];
        if (a == DP_U) Expand(energy, {st + 1, piv - 1, DP_P}, {st + 1, CTD_MISMATCH});
        Expand(energy + dp_[piv + 1][en][DP_U], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U},
            {st + 1, CTD_MISMATCH});

        // .(   ).<(   ) > Left coax - U, U2
        energy = base_energy + base11 + em_->MismatchCoaxial(pl1b, pb, stb, st1b);
        Expand(energy + dp_[piv + 1][en][DP_U_WC], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U_WC},
            {st + 1, CTD_LCOAX_WITH_NEXT}, {piv + 1, CTD_LCOAX_WITH_PREV});
        Expand(energy + dp_[piv + 1][en][DP_U_GU], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U_GU},
            {st + 1, CTD_LCOAX_WITH_NEXT}, {piv + 1, CTD_LCOAX_WITH_PREV});

        // (   )(<   ) > Flush coax - U, U2
        energy =
            base_energy + base01 + em_->stack[pl1b][pb][WcPair(pb)][stb] + dp_[piv][en][DP_U_WC];
        Expand(energy, {st, piv - 1, DP_P}, {piv, en, DP_U_WC}, {st, CTD_FCOAX_WITH_NEXT},
            {piv, CTD_FCOAX_WITH_PREV});

        if (IsGu(pb)) {
          energy =
              base_energy + base01 + em_->stack[pl1b][pb][GuPair(pb)][stb] + dp_[piv][en][DP_U_GU];
          Expand(energy, {st, piv - 1, DP_P}, {piv, en, DP_U_GU}, {st, CTD_FCOAX_WITH_NEXT},
              {piv, CTD_FCOAX_WITH_PREV});
        }

        if (piv < en - 1) {
          // (   )<.(   ). > Right coax forward - U, U2
          energy = base_energy + base00 + dp_[piv + 1][en][DP_U_RC];
          Expand(energy, {st, piv, DP_P}, {piv + 1, en, DP_U_RC}, {st, CTD_RC_WITH_NEXT},
              {piv + 2, CTD_RC_WITH_PREV});
        }
      }
    }
  }
  for (const auto& struc : finished_) {
    assert(struc.not_yet_expanded.empty());
    fn(struc.res);
  }
  return static_cast<int>(finished_.size());
}

}  // namespace mrna::subopt
