// Copyright 2016 E.
#include "fold/suboptimal1.h"

#include <algorithm>
#include <utility>
#include <vector>

#include "energy/energy_globals.h"

namespace mrna {
namespace fold {
namespace internal {

using energy::gem;
using energy::gpc;

int Suboptimal1::Run(SuboptimalCallback fn, bool sorted) {
  memset(gp.data(), -1, gp.size());
  memset(gctd.data(), CTD_NA, gctd.size());
  q.reserve(gr.size());  // Reasonable reservation.
  cache.Reserve(gr.size());

  // If require sorted output, or limited number of structures (requires sorting).
  if (sorted || max_structures != MAX_STRUCTURES) {
    int num_structures = 0;
    energy_t cur_delta = 0;
    while (num_structures < max_structures && cur_delta != MAX_E && cur_delta <= delta) {
      auto res = RunInternal(fn, cur_delta, true, max_structures - num_structures);
      num_structures += res.first;
      cur_delta = res.second;
    }
    return num_structures;
  }
  return RunInternal(fn, delta, false, MAX_STRUCTURES).first;
}

std::pair<int, int> Suboptimal1::RunInternal(
    SuboptimalCallback fn, energy_t cur_delta, bool exact_energy, int structure_limit) {
  // General idea is perform a dfs of the expand tree. Keep track of the current partial structures
  // and energy. Also keep track of what is yet to be expanded. Each node is either a terminal,
  // or leads to one expansion (either from unexpanded, or from expanding itself) - if there is
  // another expansion it is put on the unexpanded list. Everything in unexpanded never affects
  // the CTDs or energy of the current state - that is rolled into the modification when that
  // unexpanded is originally generated.

  int num_structures = 0;
  // Store the smallest energy above cur_delta we see. If we reach our |structure_limit| before
  // finishing, we might not see the smallest one, but it's okay since we won't be called again.
  // Otherwise, we will completely finish, and definitely see it.
  energy_t next_seen = MAX_E;
  energy_t energy = 0;
  q.clear();
  unexpanded.clear();
  grep.resize(gr.size(), '.');
  q.push_back({0, {0, -1, EXT}, false});
  while (!q.empty()) {
    auto& s = q.back();
    assert(s.expand.st != -1);

    const auto& exps = GetExpansion(s.expand);
    assert(!exps.empty());  // Must produce at least one expansion: {-1, -1, -1}.

    // Undo previous child's ctds and energy. The pairing is undone by the child.
    // Also remove from unexpanded if the previous child added stuff to it.
    if (s.idx != 0) {
      const auto& pexp = exps[s.idx - 1];
      if (pexp.ctd0.idx != -1) gctd[pexp.ctd0.idx] = CTD_NA;
      if (pexp.ctd1.idx != -1) gctd[pexp.ctd1.idx] = CTD_NA;
      if (pexp.unexpanded.st != -1) unexpanded.pop_back();
      energy -= pexp.energy;
    }

    // Update the next best seen variable
    if (s.idx != static_cast<int>(exps.size()) && exps[s.idx].energy + energy > cur_delta)
      next_seen = std::min(next_seen, exps[s.idx].energy + energy);

    // If we ran out of expansions, or the next expansion would take us over the delta limit
    // we are done with this node.
    if (s.idx == static_cast<int>(exps.size()) || exps[s.idx].energy + energy > cur_delta) {
      // Finished looking at this node, so undo this node's modifications to the global state.
      if (s.expand.en != -1 && s.expand.a == DP_P) {
        gp[s.expand.st] = gp[s.expand.en] = -1;
        grep[s.expand.st] = grep[s.expand.en] = '.';
      }
      if (s.should_unexpand) unexpanded.push_back(s.expand);
      q.pop_back();
      continue;  // Done.
    }

    const auto& exp = exps[s.idx++];
    dfs_state_t ns = {0, exp.to_expand, false};
    energy += exp.energy;
    if (exp.to_expand.st == -1) {
      // Can't have an unexpanded without a to_expand. Also can't set ctds or affect energy.
      assert(exp.unexpanded.st == -1);
      assert(exp.ctd0.idx == -1 && exp.ctd1.idx == -1);
      // Use an unexpanded now, if one exists.
      if (unexpanded.empty()) {
        // At a terminal state.
        if (!exact_energy || energy == cur_delta) {
          computed_t tmp_computed = {
              {std::move(gr), std::move(gp)}, std::move(gctd), energy + gext[0][EXT]};
          fn(tmp_computed);
          ++num_structures;
          // Move everything back
          gr = std::move(tmp_computed.s.r);
          gp = std::move(tmp_computed.s.p);
          gctd = std::move(tmp_computed.base_ctds);

          // Hit structure limit.
          if (num_structures == structure_limit) return {num_structures, -1};
        }
        continue;  // Done
      } else {
        ns.expand = unexpanded.back();
        unexpanded.pop_back();
        // This node should replace itself into |unexpanded| when its done.
        ns.should_unexpand = true;
      }
    } else {
      // Apply child's modifications to the global state.
      if (exp.ctd0.idx != -1) gctd[exp.ctd0.idx] = exp.ctd0.ctd;
      if (exp.ctd1.idx != -1) gctd[exp.ctd1.idx] = exp.ctd1.ctd;
      if (exp.unexpanded.st != -1) unexpanded.push_back(exp.unexpanded);
    }
    if (ns.expand.en != -1 && ns.expand.a == DP_P) {
      gp[ns.expand.st] = ns.expand.en;
      gp[ns.expand.en] = ns.expand.st;
      grep[ns.expand.st] = '(';
      grep[ns.expand.en] = ')';
    }
    q.push_back(ns);
  }
  assert(unexpanded.empty() && energy == 0 && gp == std::vector<int>(gp.size(), -1) &&
      gctd == std::vector<Ctd>(gctd.size(), CTD_NA));
  return {num_structures, next_seen};
}

std::vector<expand_t> GenerateExpansions(const index_t& to_expand, energy_t delta) {
  const int N = static_cast<int>(gr.size());
  int st = to_expand.st, en = to_expand.en, a = to_expand.a;
  std::vector<expand_t> exps;
  // Temporary variable to hold energy calculations.
  energy_t energy = 0;
  // Exterior loop
  if (en == -1) {
    if (a == EXT) {
      // Base case: do nothing.
      if (st == N)
        exps.emplace_back(0);
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

      // (   )<.( * ). > Right coax backward
      if (a == EXT_RCOAX) {
        energy = base11 + gem.MismatchCoaxial(en1b, enb, stb, st1b) + gext[en + 1][EXT];
        // We don't set ctds here, since we already set them in the forward case.
        if (energy <= delta) exps.push_back({energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P}});
      }

      // Cases for EXT, EXT_WC, EXT_GU.
      // (   )<   >
      // If we are at EXT then this is unused.
      energy = base00 + gext[en + 1][EXT];
      if (energy <= delta) {
        if (a == EXT) exps.push_back({energy, {en + 1, -1, EXT}, {st, en, DP_P}, {st, CTD_UNUSED}});

        // (   )<   >
        // If we are at EXT_WC or EXT_GU, the CTDs for this have already have been set from a
        // coaxial stack.
        if ((a == EXT_WC && IsWatsonCrick(stb, enb)) || (a == EXT_GU && IsGu(stb, enb)))
          exps.push_back({energy, {en + 1, -1, EXT}, {st, en, DP_P}});
      }

      // Everything after this is only for EXT.
      if (a != EXT) continue;

      // (   )3<   > 3'
      energy = base01 + gem.dangle3[en1b][enb][stb] + gext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({energy, {en + 1, -1, EXT}, {st, en - 1, DP_P}, {st, CTD_3_DANGLE}});

      // 5(   )<   > 5'
      energy = base10 + gem.dangle5[enb][stb][st1b] + gext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({energy, {en + 1, -1, EXT}, {st + 1, en, DP_P}, {st + 1, CTD_5_DANGLE}});

      // .(   ).<   > Terminal mismatch
      energy = base11 + gem.terminal[en1b][enb][stb][st1b] + gext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P}, {st + 1, CTD_MISMATCH}});

      if (en < N - 1) {
        // .(   ).<(   ) > Left coax
        energy = base11 + gem.MismatchCoaxial(en1b, enb, stb, st1b);
        if (energy + gext[en + 1][EXT_GU] <= delta)
          exps.push_back(
              {energy + gext[en + 1][EXT_GU], {en + 1, -1, EXT_GU}, {st + 1, en - 1, DP_P},
                  {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT}});
        if (energy + gext[en + 1][EXT_WC] <= delta)
          exps.push_back(
              {energy + gext[en + 1][EXT_WC], {en + 1, -1, EXT_WC}, {st + 1, en - 1, DP_P},
                  {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT}});
      }

      if (en < N - 2) {
        // (   )<.(   ). > Right coax forward
        energy = base00 + gext[en + 1][EXT_RCOAX];
        if (energy <= delta)
          exps.push_back({energy, {en + 1, -1, EXT_RCOAX}, {st, en, DP_P},
              {en + 2, CTD_RCOAX_WITH_PREV}, {st, CTD_RCOAX_WITH_NEXT}});
      }

      // (   )(<   ) > Flush coax
      energy = base01 + gem.stack[en1b][enb][enb ^ 3][stb] + gext[en][EXT_WC];
      if (energy <= delta)
        exps.push_back({energy, {en, -1, EXT_WC}, {st, en - 1, DP_P}, {en, CTD_FCOAX_WITH_PREV},
            {st, CTD_FCOAX_WITH_NEXT}});

      if (enb == G || enb == U) {
        energy = base01 + gem.stack[en1b][enb][enb ^ 1][stb] + gext[en][EXT_GU];
        if (energy <= delta)
          exps.push_back({energy, {en, -1, EXT_GU}, {st, en - 1, DP_P}, {en, CTD_FCOAX_WITH_PREV},
              {st, CTD_FCOAX_WITH_NEXT}});
      }
    }
    // Finished exterior loop, don't do anymore.
    return exps;
  }

  // Declare the usual base aliases.
  const auto stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en], en1b = gr[en - 1],
             en2b = gr[en - 2];

  // Normal stuff
  if (a == DP_P) {
    // Two loops.
    int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
    for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
      for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
        energy = energy::FastTwoLoop(st, en, ist, ien) + gdp[ist][ien][DP_P] - gdp[st][en][a];
        if (energy <= delta) exps.push_back({energy, {ist, ien, DP_P}});
      }
    }

    // Hairpin loop
    energy = energy::FastHairpin(st, en) - gdp[st][en][a];
    if (energy <= delta) exps.emplace_back(energy);

    auto base_and_branch = gpc.augubranch[stb][enb] + gem.multiloop_hack_a - gdp[st][en][a];
    // (<   ><    >)
    energy = base_and_branch + gdp[st + 1][en - 1][DP_U2];
    if (energy <= delta) exps.push_back({energy, {st + 1, en - 1, DP_U2}, {en, CTD_UNUSED}});
    // (3<   ><   >) 3'
    energy = base_and_branch + gdp[st + 2][en - 1][DP_U2] + gem.dangle3[stb][st1b][enb];
    if (energy <= delta) exps.push_back({energy, {st + 2, en - 1, DP_U2}, {en, CTD_3_DANGLE}});
    // (<   ><   >5) 5'
    energy = base_and_branch + gdp[st + 1][en - 2][DP_U2] + gem.dangle5[stb][en1b][enb];
    if (energy <= delta) exps.push_back({energy, {st + 1, en - 2, DP_U2}, {en, CTD_5_DANGLE}});
    // (.<   ><   >.) Terminal mismatch
    energy = base_and_branch + gdp[st + 2][en - 2][DP_U2] + gem.terminal[stb][st1b][en1b][enb];
    if (energy <= delta) exps.push_back({energy, {st + 2, en - 2, DP_U2}, {en, CTD_MISMATCH}});

    for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
      base_t pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1], pr1b = gr[piv + 2];

      // (.(   )   .) Left outer coax - P
      auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);
      energy = base_and_branch + gdp[st + 2][piv][DP_P] + gpc.augubranch[st2b][plb] +
          gdp[piv + 1][en - 2][DP_U] + outer_coax;
      if (energy <= delta)
        exps.push_back({energy, {st + 2, piv, DP_P}, {piv + 1, en - 2, DP_U},
            {st + 2, CTD_LCOAX_WITH_PREV}, {en, CTD_LCOAX_WITH_NEXT}});

      // (.   (   ).) Right outer coax
      energy = base_and_branch + gdp[st + 2][piv][DP_U] + gpc.augubranch[prb][en2b] +
          gdp[piv + 1][en - 2][DP_P] + outer_coax;
      if (energy <= delta)
        exps.push_back({energy, {st + 2, piv, DP_U}, {piv + 1, en - 2, DP_P},
            {piv + 1, CTD_RCOAX_WITH_NEXT}, {en, CTD_RCOAX_WITH_PREV}});

      // (.(   ).   ) Left right coax
      energy = base_and_branch + gdp[st + 2][piv - 1][DP_P] + gpc.augubranch[st2b][pl1b] +
          gdp[piv + 1][en - 1][DP_U] + gem.MismatchCoaxial(pl1b, plb, st1b, st2b);
      if (energy <= delta)
        exps.push_back({energy, {st + 2, piv - 1, DP_P}, {piv + 1, en - 1, DP_U},
            {st + 2, CTD_RCOAX_WITH_PREV}, {en, CTD_RCOAX_WITH_NEXT}});

      // (   .(   ).) Right left coax
      energy = base_and_branch + gdp[st + 1][piv][DP_U] + gpc.augubranch[pr1b][en2b] +
          gdp[piv + 2][en - 2][DP_P] + gem.MismatchCoaxial(en2b, en1b, prb, pr1b);
      if (energy <= delta)
        exps.push_back({energy, {st + 1, piv, DP_U}, {piv + 2, en - 2, DP_P},
            {piv + 2, CTD_LCOAX_WITH_NEXT}, {en, CTD_LCOAX_WITH_PREV}});

      // ((   )   ) Left flush coax
      energy = base_and_branch + gdp[st + 1][piv][DP_P] + gpc.augubranch[st1b][plb] +
          gdp[piv + 1][en - 1][DP_U] + gem.stack[stb][st1b][plb][enb];
      if (energy <= delta)
        exps.push_back({energy, {st + 1, piv, DP_P}, {piv + 1, en - 1, DP_U},
            {st + 1, CTD_FCOAX_WITH_PREV}, {en, CTD_FCOAX_WITH_NEXT}});

      // (   (   )) Right flush coax
      energy = base_and_branch + gdp[st + 1][piv][DP_U] + gpc.augubranch[prb][en1b] +
          gdp[piv + 1][en - 1][DP_P] + gem.stack[stb][prb][en1b][enb];
      if (energy <= delta)
        exps.push_back({energy, {st + 1, piv, DP_U}, {piv + 1, en - 1, DP_P},
            {piv + 1, CTD_FCOAX_WITH_NEXT}, {en, CTD_FCOAX_WITH_PREV}});
    }
    return exps;
  }

  // Left unpaired. Either DP_U or DP_U2.
  if (st + 1 < en && (a == DP_U || a == DP_U2)) {
    energy = gdp[st + 1][en][a] - gdp[st][en][a];
    if (energy <= delta) exps.push_back({energy, {st + 1, en, a}});
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
    // (   )<.( ** ). > Right coax backward
    if (a == DP_U_RCOAX) {
      energy = base11 + gem.MismatchCoaxial(pl1b, pb, stb, st1b);
      // Our ctds will have already been set by now.
      if (energy <= delta) exps.push_back({energy, {st + 1, piv - 1, DP_P}});
      if (energy + gdp[piv + 1][en][DP_U] <= delta)
        exps.push_back(
            {energy + gdp[piv + 1][en][DP_U], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U}});
      continue;
    }
    // From here on, a must be U, U2, U_WC, or U_GU.

    // (   )<   > - U, U2, U_WC?, U_GU?
    energy = base00;
    if (a == DP_U) {
      if (energy <= delta) exps.push_back({energy, {st, piv, DP_P}, {st, CTD_UNUSED}});
      if (energy + gdp[piv + 1][en][DP_U] <= delta)
        exps.push_back({energy + gdp[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U},
            {st, CTD_UNUSED}});
    }
    if (a == DP_U2 && energy + gdp[piv + 1][en][DP_U] <= delta)
      exps.push_back({energy + gdp[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U},
          {st, CTD_UNUSED}});
    if (a == DP_U_WC || a == DP_U_GU) {
      // Make sure we don't form any branches that are not the right type of pair.
      if ((a == DP_U_WC && IsWatsonCrick(stb, pb)) || (a == DP_U_GU && IsGu(stb, pb))) {
        if (energy <= delta) exps.push_back({energy, {st, piv, DP_P}});
        if (energy + gdp[piv + 1][en][DP_U] <= delta)
          exps.push_back({energy + gdp[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U}});
      }
      continue;
    }
    // From here on, a must be U or U2.
    assert(a == DP_U || a == DP_U2);

    // (   )3<   > 3' - U, U2
    energy = base01 + gem.dangle3[pl1b][pb][stb];
    // Can only let the rest be unpaired if we only need one branch, i.e. DP_U not DP_U2.
    if (a == DP_U && energy <= delta)
      exps.push_back({energy, {st, piv - 1, DP_P}, {st, CTD_3_DANGLE}});
    if (energy + gdp[piv + 1][en][DP_U] <= delta)
      exps.push_back({energy + gdp[piv + 1][en][DP_U], {st, piv - 1, DP_P}, {piv + 1, en, DP_U},
          {st, CTD_3_DANGLE}});

    // 5(   )<   > 5' - U, U2
    energy = base10 + gem.dangle5[pb][stb][st1b];
    if (a == DP_U && energy <= delta)
      exps.push_back({energy, {st + 1, piv, DP_P}, {st + 1, CTD_5_DANGLE}});
    if (energy + gdp[piv + 1][en][DP_U] <= delta)
      exps.push_back({energy + gdp[piv + 1][en][DP_U], {st + 1, piv, DP_P}, {piv + 1, en, DP_U},
          {st + 1, CTD_5_DANGLE}});

    // .(   ).<   > Terminal mismatch - U, U2
    energy = base11 + gem.terminal[pl1b][pb][stb][st1b];
    if (a == DP_U && energy <= delta)
      exps.push_back({energy, {st + 1, piv - 1, DP_P}, {}, {st + 1, CTD_MISMATCH}});
    if (energy + gdp[piv + 1][en][DP_U] <= delta)
      exps.push_back({energy + gdp[piv + 1][en][DP_U], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U},
          {st + 1, CTD_MISMATCH}});

    // .(   ).<(   ) > Left coax - U, U2
    energy = base11 + gem.MismatchCoaxial(pl1b, pb, stb, st1b);
    if (energy + gdp[piv + 1][en][DP_U_WC] <= delta)
      exps.push_back({energy + gdp[piv + 1][en][DP_U_WC], {st + 1, piv - 1, DP_P},
          {piv + 1, en, DP_U_WC}, {st + 1, CTD_LCOAX_WITH_NEXT}, {piv + 1, CTD_LCOAX_WITH_PREV}});
    if (energy + gdp[piv + 1][en][DP_U_GU] <= delta)
      exps.push_back({energy + gdp[piv + 1][en][DP_U_GU], {st + 1, piv - 1, DP_P},
          {piv + 1, en, DP_U_GU}, {st + 1, CTD_LCOAX_WITH_NEXT}, {piv + 1, CTD_LCOAX_WITH_PREV}});

    // (   )<.(   ). > Right coax forward - U, U2
    energy = base00 + gdp[piv + 1][en][DP_U_RCOAX];
    if (energy <= delta)
      exps.push_back({energy, {st, piv, DP_P}, {piv + 1, en, DP_U_RCOAX}, {st, CTD_RCOAX_WITH_NEXT},
          {piv + 2, CTD_RCOAX_WITH_PREV}});

    // (   )(<   ) > Flush coax - U, U2
    energy = base01 + gem.stack[pl1b][pb][pb ^ 3][stb] + gdp[piv][en][DP_U_WC];
    if (energy <= delta)
      exps.push_back({energy, {st, piv - 1, DP_P}, {piv, en, DP_U_WC}, {st, CTD_FCOAX_WITH_NEXT},
          {piv, CTD_FCOAX_WITH_PREV}});

    if (pb == G || pb == U) {
      energy = base01 + gem.stack[pl1b][pb][pb ^ 1][stb] + gdp[piv][en][DP_U_GU];
      if (energy <= delta)
        exps.push_back({energy, {st, piv - 1, DP_P}, {piv, en, DP_U_GU}, {st, CTD_FCOAX_WITH_NEXT},
            {piv, CTD_FCOAX_WITH_PREV}});
    }
  }

  return exps;
}

}  // namespace internal
}  // namespace fold
}  // namespace mrna
