#include "base.h"
#include "fold/fold_internal.h"
#include "fold/suboptimal1_base.h"

namespace memerna {
namespace fold {
namespace internal {

using namespace energy;

std::vector<expand_t> GenerateExpansions(const index_t& to_expand) {
  const int N = int(gr.size());
  int st = to_expand.st, en = to_expand.en, a = to_expand.a;
  std::vector<expand_t> exps;
  // Temporary variable to hold energy calculations.
  energy_t energy = 0;
  // Exterior loop
  if (en == -1) {
    if (a == EXT) {
      // Base case: do nothing.
      if (st == N)
        exps.push_back({0});
      else
        // Case: No pair starting here (for EXT only)
        exps.push_back({gext[st + 1][EXT] - gext[st][a], {st + 1, -1, EXT}});
    }
    for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = gr[st], st1b = gr[st + 1], enb = gr[en], en1b = gr[en - 1];
      const auto base00 = gdp[st][en][DP_P] + gem.AuGuPenalty(stb, enb) - gext[st][a];
      const auto base01 = gdp[st][en - 1][DP_P] + gem.AuGuPenalty(stb, en1b) - gext[st][a];
      const auto base10 = gdp[st + 1][en][DP_P] + gem.AuGuPenalty(st1b, enb) - gext[st][a];
      const auto base11 = gdp[st + 1][en - 1][DP_P] + gem.AuGuPenalty(st1b, en1b) - gext[st][a];

      // (   ).<( * ). > Right coax backward
      if (st > 0 && a == EXT_RCOAX) {
        energy = base01 + gem.MismatchCoaxial(en1b, enb, gr[st - 1], stb) + gext[en + 1][EXT];
        // We don't set ctds here, since we already set them in the forward case.
        exps.push_back({energy, {en + 1, -1, EXT}, {st, en - 1, DP_P}});
      }

      // Cases for EXT, EXT_WC, EXT_GU.
      // (   )<   >
      // If we are at EXT then this is unused.
      energy = base00 + gext[en + 1][EXT];
      if (a == EXT)
        exps.push_back({energy, {en + 1, -1, EXT}, {st, en, DP_P}, {st, CTD_UNUSED}});

      // (   )<   >
      // If we are at EXT_WC or EXT_GU, the CTDs for this have already have been set from a coaxial stack.
      if ((a == EXT_WC && IsWatsonCrick(stb, enb)) || (a == EXT_GU && IsGu(stb, enb)))
        exps.push_back({energy, {en + 1, -1, EXT}, {st, en, DP_P}});

      // Everything after this is only for EXT.
      if (a != EXT) continue;

      // (   )3<   > 3'
      energy = base01 + gem.dangle3[en1b][enb][stb] + gext[en + 1][EXT];
      exps.push_back({energy, {en + 1, -1, EXT}, {st, en - 1, DP_P}, {st, CTD_3_DANGLE}});

      // 5(   )<   > 5'
      energy = base10 + gem.dangle5[enb][stb][st1b] + gext[en + 1][EXT];
      exps.push_back({energy, {en + 1, -1, EXT}, {st + 1, en, DP_P}, {st + 1, CTD_5_DANGLE}});

      // .(   ).<   > Terminal mismatch
      energy = base11 + gem.terminal[en1b][enb][stb][st1b] + gext[en + 1][EXT];
      exps.push_back({energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P},
          {st + 1, CTD_MISMATCH}});

      if (en < N - 1) {
        // .(   ).<(   ) > Left coax
        energy = base11 + gem.MismatchCoaxial(en1b, enb, stb, st1b);
        exps.push_back({energy + gext[en + 1][EXT_GU], {en + 1, -1, EXT_GU}, {st + 1, en - 1, DP_P},
            {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT}});
        exps.push_back({energy + gext[en + 1][EXT_WC], {en + 1, -1, EXT_WC}, {st + 1, en - 1, DP_P},
            {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT}});

        // (   ).<(   ). > Right coax forward
        energy = base01 + gext[en + 1][EXT_RCOAX];
        exps.push_back({energy, {en + 1, -1, EXT_RCOAX}, {st, en - 1, DP_P},
            {en + 1, CTD_RCOAX_WITH_PREV}, {st, CTD_RCOAX_WITH_NEXT}});

        // (   )<(   ) > Flush coax
        const auto enrb = gr[en + 1];
        energy = base00 + gem.stack[enb][enrb][enrb ^ 3][stb] + gext[en + 1][EXT_WC];
        exps.push_back({energy, {en + 1, -1, EXT_WC}, {st, en, DP_P},
            {en + 1, CTD_FCOAX_WITH_PREV}, {st, CTD_FCOAX_WITH_NEXT}});

        if (enrb == G || enrb == U) {
          energy = base00 + gem.stack[enb][enrb][enrb ^ 1][stb] + gext[en + 1][EXT_GU];
          exps.push_back({energy, {en + 1, -1, EXT_GU}, {st, en, DP_P},
              {en + 1, CTD_FCOAX_WITH_PREV}, {st, CTD_FCOAX_WITH_NEXT}});
        }
      }
    }
    // Finished exterior loop, don't do anymore.
    return exps;
  }

  // Declare the usual base aliases.
  const auto stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en], en1b = gr[en - 1], en2b = gr[en - 2];

  // Normal stuff
  if (a == DP_P) {
    // Two loops.
    int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
    for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
      for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
        energy = FastTwoLoop(st, en, ist, ien) + gdp[ist][ien][DP_P] - gdp[st][en][a];
        exps.push_back({energy, {ist, ien, DP_P}});
      }
    }

    // Hairpin loop
    energy = FastHairpin(st, en) - gdp[st][en][a];
    exps.push_back({energy});

    auto base_and_branch = gpc.augubranch[stb][enb] + gem.multiloop_hack_a - gdp[st][en][a];
    // (<   ><    >)
    energy = base_and_branch + gdp[st + 1][en - 1][DP_U2];
    exps.push_back({energy, {st + 1, en - 1, DP_U2}, {en, CTD_UNUSED}});
    // (3<   ><   >) 3'
    energy = base_and_branch + gdp[st + 2][en - 1][DP_U2] + gem.dangle3[stb][st1b][enb];
    exps.push_back({energy, {st + 2, en - 1, DP_U2}, {en, CTD_3_DANGLE}});
    // (<   ><   >5) 5'
    energy = base_and_branch + gdp[st + 1][en - 2][DP_U2] + gem.dangle5[stb][en1b][enb];
    exps.push_back({energy, {st + 1, en - 2, DP_U2}, {en, CTD_5_DANGLE}});
    // (.<   ><   >.) Terminal mismatch
    energy = base_and_branch + gdp[st + 2][en - 2][DP_U2] + gem.terminal[stb][st1b][en1b][enb];
    exps.push_back({energy, {st + 2, en - 2, DP_U2}, {en, CTD_MISMATCH}});

    for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
      base_t pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1], pr1b = gr[piv + 2];

      // (.(   )   .) Left outer coax - P
      auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);
      energy = base_and_branch + gdp[st + 2][piv][DP_P] + gpc.augubranch[st2b][plb] +
          gdp[piv + 1][en - 2][DP_U] + outer_coax;
      exps.push_back({energy, {st + 2, piv, DP_P}, {piv + 1, en - 2, DP_U},
          {st + 2, CTD_LCOAX_WITH_PREV}, {en, CTD_LCOAX_WITH_NEXT}});

      // (.   (   ).) Right outer coax
      energy = base_and_branch + gdp[st + 2][piv][DP_U] + gpc.augubranch[prb][en2b] +
          gdp[piv + 1][en - 2][DP_P] + outer_coax;
      exps.push_back({energy, {st + 2, piv, DP_U}, {piv + 1, en - 2, DP_P},
          {piv + 1, CTD_RCOAX_WITH_NEXT}, {en, CTD_RCOAX_WITH_PREV}});

      // (.(   ).   ) Left right coax
      energy = base_and_branch + gdp[st + 2][piv - 1][DP_P] + gpc.augubranch[st2b][pl1b] +
          gdp[piv + 1][en - 1][DP_U] + gem.MismatchCoaxial(pl1b, plb, st1b, st2b);
      exps.push_back({energy, {st + 2, piv - 1, DP_P}, {piv + 1, en - 1, DP_U},
          {st + 2, CTD_RCOAX_WITH_PREV}, {en, CTD_RCOAX_WITH_NEXT}});

      // (   .(   ).) Right left coax
      energy = base_and_branch + gdp[st + 1][piv][DP_U] + gpc.augubranch[pr1b][en2b] +
          gdp[piv + 2][en - 2][DP_P] + gem.MismatchCoaxial(en2b, en1b, prb, pr1b);
      exps.push_back({energy, {st + 1, piv, DP_U}, {piv + 2, en - 2, DP_P},
          {piv + 2, CTD_LCOAX_WITH_NEXT}, {en, CTD_LCOAX_WITH_PREV}});

      // ((   )   ) Left flush coax
      energy = base_and_branch + gdp[st + 1][piv][DP_P] + gpc.augubranch[st1b][plb] +
          gdp[piv + 1][en - 1][DP_U] + gem.stack[stb][st1b][plb][enb];
      exps.push_back({energy, {st + 1, piv, DP_P}, {piv + 1, en - 1, DP_U},
          {st + 1, CTD_FCOAX_WITH_PREV}, {en, CTD_FCOAX_WITH_NEXT}});

      // (   (   )) Right flush coax
      energy = base_and_branch + gdp[st + 1][piv][DP_U] + gpc.augubranch[prb][en1b] +
          gdp[piv + 1][en - 1][DP_P] + gem.stack[stb][prb][en1b][enb];
      exps.push_back({energy, {st + 1, piv, DP_U}, {piv + 1, en - 1, DP_P},
          {piv + 1, CTD_FCOAX_WITH_NEXT}, {en, CTD_FCOAX_WITH_PREV}});
    }
    return exps;
  }

  // Left unpaired. Either DP_U or DP_U2.
  if (st + 1 < en && (a == DP_U || a == DP_U2)) {
    energy = gdp[st + 1][en][a] - gdp[st][en][a];
    exps.push_back({energy, {st + 1, en, a}});
  }

  // Pair here.
  for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
    //   (   .   )<   (
    // stb pl1b pb   pr1b
    auto pb = gr[piv], pl1b = gr[piv - 1];
    // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
    auto base00 = gdp[st][piv][DP_P] + gpc.augubranch[stb][pb] - gdp[st][en][a];
    auto base01 = gdp[st][piv - 1][DP_P] + gpc.augubranch[stb][pl1b] - gdp[st][en][a];
    auto base10 = gdp[st + 1][piv][DP_P] + gpc.augubranch[st1b][pb] - gdp[st][en][a];
    auto base11 = gdp[st + 1][piv - 1][DP_P] + gpc.augubranch[st1b][pl1b] - gdp[st][en][a];

    // Check a == U_RCOAX:
    // (   ).<( ** ). > Right coax backward
    if (a == DP_U_RCOAX) {
      if (st > 0) {
        energy = base01 + gem.MismatchCoaxial(pl1b, pb, gr[st - 1], stb);
        // Our ctds will have already been set by now.
        exps.push_back({energy, {st, piv - 1, DP_P}});
        exps.push_back({energy + gdp[piv + 1][en][DP_U], {st, piv - 1, DP_P},
            {piv + 1, en, DP_U}});
      }
      continue;
    }
    // From here on, a must be U, U2, U_WC, or U_GU.

    // (   )<   > - U, U2, U_WC?, U_GU?
    energy = base00;
    if (a == DP_U) {
      exps.push_back({energy, {st, piv, DP_P}, {st, CTD_UNUSED}});
      exps.push_back({energy + gdp[piv + 1][en][DP_U], {st, piv, DP_P},
          {piv + 1, en, DP_U}, {st, CTD_UNUSED}});
    }
    if (a == DP_U2)
      exps.push_back({energy + gdp[piv + 1][en][DP_U], {st, piv, DP_P},
          {piv + 1, en, DP_U}, {st, CTD_UNUSED}});
    if (a == DP_U_WC || a == DP_U_GU) {
      // Make sure we don't form any branches that are not the right type of pair.
      if ((a == DP_U_WC && IsWatsonCrick(stb, pb)) || (a == DP_U_GU && IsGu(stb, pb))) {
        exps.push_back({energy, {st, piv, DP_P}});
        exps.push_back({energy + gdp[piv + 1][en][DP_U], {st, piv, DP_P},
            {piv + 1, en, DP_U}});
      }
      continue;
    }
    // From here on, a must be U or U2.
    assert(a == DP_U || a == DP_U2);

    // (   )3<   > 3' - U, U2
    energy = base01 + gem.dangle3[pl1b][pb][stb];
    // Can only let the rest be unpaired if we only need one branch, i.e. DP_U not DP_U2.
    if (a == DP_U) exps.push_back({energy, {st, piv - 1, DP_P}, {st, CTD_3_DANGLE}});
    exps.push_back({energy + gdp[piv + 1][en][DP_U],
        {st, piv - 1, DP_P}, {piv + 1, en, DP_U}, {st, CTD_3_DANGLE}});

    // 5(   )<   > 5' - U, U2
    energy = base10 + gem.dangle5[pb][stb][st1b];
    if (a == DP_U) exps.push_back({energy, {st + 1, piv, DP_P}, {st + 1, CTD_5_DANGLE}});
    exps.push_back({energy + gdp[piv + 1][en][DP_U],
        {st + 1, piv, DP_P}, {piv + 1, en, DP_U}, {st + 1, CTD_5_DANGLE}});

    // .(   ).<   > Terminal mismatch - U, U2
    energy = base11 + gem.terminal[pl1b][pb][stb][st1b];
    if (a == DP_U)
      exps.push_back({energy, {st + 1, piv - 1, DP_P},
          {}, {st + 1, CTD_MISMATCH}});
    exps.push_back({energy + gdp[piv + 1][en][DP_U], {st + 1, piv - 1, DP_P},
        {piv + 1, en, DP_U}, {st + 1, CTD_MISMATCH}});

    // .(   ).<(   ) > Left coax - U, U2
    energy = base11 + gem.MismatchCoaxial(pl1b, pb, stb, st1b);
    exps.push_back({energy + gdp[piv + 1][en][DP_U_WC], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U_WC},
        {st + 1, CTD_LCOAX_WITH_NEXT}, {piv + 1, CTD_LCOAX_WITH_PREV}});
    exps.push_back({energy + gdp[piv + 1][en][DP_U_GU], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U_GU},
        {st + 1, CTD_LCOAX_WITH_NEXT}, {piv + 1, CTD_LCOAX_WITH_PREV}});

    // (   ).<(   ). > Right coax forward - U, U2
    energy = base01 + gdp[piv + 1][en][DP_U_RCOAX];
    exps.push_back({energy, {st, piv - 1, DP_P}, {piv + 1, en, DP_U_RCOAX},
        {st, CTD_RCOAX_WITH_NEXT}, {piv + 1, CTD_RCOAX_WITH_PREV}});

    // There has to be remaining bases to even have a chance at these cases.
    if (piv < en) {
      auto pr1b = gr[piv + 1];
      // (   )<(   ) > Flush coax - U, U2
      energy = base00 + gem.stack[pb][pr1b][pr1b ^ 3][stb] + gdp[piv + 1][en][DP_U_WC];
      exps.push_back({energy, {st, piv, DP_P}, {piv + 1, en, DP_U_WC},
          {st, CTD_FCOAX_WITH_NEXT}, {piv + 1, CTD_FCOAX_WITH_PREV}});

      if (pr1b == G || pr1b == U) {
        energy = base00 + gem.stack[pb][pr1b][pr1b ^ 1][stb] + gdp[piv + 1][en][DP_U_GU];
        exps.push_back({energy, {st, piv, DP_P}, {piv + 1, en, DP_U_GU},
            {st, CTD_FCOAX_WITH_NEXT}, {piv + 1, CTD_FCOAX_WITH_PREV}});
      }
    }
  }

  return exps;
}

}
}
}
