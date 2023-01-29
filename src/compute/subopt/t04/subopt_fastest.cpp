// Copyright 2016 Eliot Courtney.
#include "compute/subopt/t04/subopt_fastest.h"

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "compute/subopt/subopt.h"
#include "compute/traceback/traceback.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/secondary.h"

namespace mrna::subopt::t04 {

SuboptFastest::SuboptFastest(
    Primary r, erg::t04::Model::Ptr em, DpArray dp, ExtArray ext, SuboptCfg cfg)
    : r_(std::move(r)), em_(std::move(em)), pc_(Primary(r_), em_), dp_(std::move(dp)),
      ext_(std::move(ext)), cfg_(cfg) {}

int SuboptFastest::Run(const SuboptCallback& fn) {
  res_ = SuboptResult(ZERO_E, tb::TracebackResult(Secondary(r_.size()), Ctds(r_.size())));
  q_.reserve(r_.size());  // Reasonable reservation.
  cache_.Reserve(r_.size());

  // If require sorted output, or limited number of structures (requires sorting).
  if (cfg_.sorted || cfg_.strucs != SuboptCfg::MAX_STRUCTURES) {
    int count = 0;
    Energy delta = ZERO_E;
    while (count < cfg_.strucs && delta != MAX_E && delta <= cfg_.delta) {
      auto res = RunInternal(fn, delta, true, cfg_.strucs - count);
      count += res.first;
      delta = res.second;
    }
    return count;
  }
  return RunInternal(fn, cfg_.delta, false, cfg_.strucs).first;
}

std::pair<int, Energy> SuboptFastest::RunInternal(
    const SuboptCallback& fn, Energy delta, bool exact_energy, int max) {
  // General idea is perform a dfs of the expand tree. Keep track of the current partial structures
  // and energy. Also keep track of what is yet to be expanded. Each node is either a terminal,
  // or leads to one expansion (either from unexpanded, or from expanding itself) - if there is
  // another expansion it is put on the unexpanded list. Everything in unexpanded never affects
  // the CTDs or energy of the current state - that is rolled into the modification when that
  // unexpanded is originally generated.

  int count = 0;
  // Store the smallest energy above delta we see. If we reach our |structure_limit| before
  // finishing, we might not see the smallest one, but it's okay since we won't be called again.
  // Otherwise, we will completely finish, and definitely see it.
  Energy next_seen = MAX_E;
  Energy energy = ZERO_E;
  q_.clear();
  unexpanded_.clear();
  q_.push_back({0, {0, -1, EXT}, false});
  while (!q_.empty()) {
    auto& s = q_.back();
    assert(s.expand.st != -1);

    const auto& exps = GetExpansion(s.expand);
    assert(!exps.empty());  // Must produce at least one expansion: {-1, -1, -1}.

    // Undo previous child's ctds and energy. The pairing is undone by the child.
    // Also remove from unexpanded if the previous child added stuff to it.
    if (s.idx != 0) {
      const auto& pexp = exps[s.idx - 1];
      if (pexp.ctd0.idx != -1) res_.tb.ctd[pexp.ctd0.idx] = CTD_NA;
      if (pexp.ctd1.idx != -1) res_.tb.ctd[pexp.ctd1.idx] = CTD_NA;
      if (pexp.unexpanded.st != -1) unexpanded_.pop_back();
      energy -= pexp.energy;
    }

    // Update the next best seen variable
    if (s.idx != static_cast<int>(exps.size()) && exps[s.idx].energy + energy > delta)
      next_seen = std::min(next_seen, exps[s.idx].energy + energy);

    // If we ran out of expansions, or the next expansion would take us over the delta limit
    // we are done with this node.
    if (s.idx == static_cast<int>(exps.size()) || exps[s.idx].energy + energy > delta) {
      // Finished looking at this node, so undo this node's modifications to the global state.
      if (s.expand.en != -1 && s.expand.a == DP_P)
        res_.tb.s[s.expand.st] = res_.tb.s[s.expand.en] = -1;
      if (s.should_unexpand) unexpanded_.push_back(s.expand);
      q_.pop_back();
      continue;  // Done.
    }

    const auto& ex = exps[s.idx++];
    DfsState ns = {0, ex.to_expand, false};
    energy += ex.energy;
    if (ex.to_expand.st == -1) {
      // Can't have an unexpanded without a to_expand. Also can't set ctds or affect energy.
      assert(ex.unexpanded.st == -1);
      assert(ex.ctd0.idx == -1 && ex.ctd1.idx == -1);
      // Use an unexpanded now, if one exists.
      if (unexpanded_.empty()) {
        // At a terminal state.
        if (!exact_energy || energy == delta) {
          res_.energy = energy + ext_[0][EXT];
          fn(res_);
          ++count;

          // Hit structure limit.
          if (count == max) return {count, CAP_E};
        }
        continue;  // Done
      }
      ns.expand = unexpanded_.back();
      unexpanded_.pop_back();
      // This node should replace itself into |unexpanded| when its done.
      ns.should_unexpand = true;

    } else {
      // Apply child's modifications to the global state.
      if (ex.ctd0.idx != -1) res_.tb.ctd[ex.ctd0.idx] = ex.ctd0.ctd;
      if (ex.ctd1.idx != -1) res_.tb.ctd[ex.ctd1.idx] = ex.ctd1.ctd;
      if (ex.unexpanded.st != -1) unexpanded_.push_back(ex.unexpanded);
    }
    if (ns.expand.en != -1 && ns.expand.a == DP_P) {
      res_.tb.s[ns.expand.st] = ns.expand.en;
      res_.tb.s[ns.expand.en] = ns.expand.st;
    }
    q_.push_back(ns);
  }
  assert(unexpanded_.empty() && energy == ZERO_E && res_.tb.s == Secondary(res_.tb.s.size()) &&
      res_.tb.ctd == Ctds(res_.tb.ctd.size()));
  return {count, next_seen};
}

std::vector<Expand> SuboptFastest::GenerateExpansions(const Index& to_expand, Energy delta) const {
  const int N = static_cast<int>(r_.size());
  const int st = to_expand.st;
  int en = to_expand.en;
  const int a = to_expand.a;
  std::vector<Expand> exps;
  // Temporary variable to hold energy calculations.
  Energy energy = ZERO_E;
  // Exterior loop
  if (en == -1) {
    if (a == EXT) {
      // Base case: do nothing.
      if (st == N)
        exps.emplace_back(ZERO_E);
      else
        // Case: No pair starting here (for EXT only)
        exps.push_back({ext_[st + 1][EXT] - ext_[st][a], {st + 1, -1, EXT}});
    }
    for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = r_[st];
      const auto st1b = r_[st + 1];
      const auto enb = r_[en];
      const auto en1b = r_[en - 1];
      const auto base00 = dp_[st][en][DP_P] + em_->AuGuPenalty(stb, enb) - ext_[st][a];
      const auto base01 = dp_[st][en - 1][DP_P] + em_->AuGuPenalty(stb, en1b) - ext_[st][a];
      const auto base10 = dp_[st + 1][en][DP_P] + em_->AuGuPenalty(st1b, enb) - ext_[st][a];
      const auto base11 = dp_[st + 1][en - 1][DP_P] + em_->AuGuPenalty(st1b, en1b) - ext_[st][a];

      // (   )<.( * ). > Right coax backward
      if (a == EXT_RC) {
        energy = base11 + em_->MismatchCoaxial(en1b, enb, stb, st1b) + ext_[en + 1][EXT];
        // We don't set ctds here, since we already set them in the forward case.
        if (energy <= delta) exps.push_back({energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P}});
      }

      // Cases for EXT, EXT_WC, EXT_GU.
      // (   )<   >
      // If we are at EXT then this is unused.
      energy = base00 + ext_[en + 1][EXT];
      if (energy <= delta) {
        if (a == EXT) exps.push_back({energy, {en + 1, -1, EXT}, {st, en, DP_P}, {st, CTD_UNUSED}});

        // (   )<   >
        // If we are at EXT_WC or EXT_GU, the CTDs for this have already have been set from a
        // coaxial stack.
        if ((a == EXT_WC && IsWcPair(stb, enb)) || (a == EXT_GU && IsGuPair(stb, enb)))
          exps.push_back({energy, {en + 1, -1, EXT}, {st, en, DP_P}});
      }

      // Everything after this is only for EXT.
      if (a != EXT) continue;

      // (   )3<   > 3'
      energy = base01 + em_->dangle3[en1b][enb][stb] + ext_[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({energy, {en + 1, -1, EXT}, {st, en - 1, DP_P}, {st, CTD_3_DANGLE}});

      // 5(   )<   > 5'
      energy = base10 + em_->dangle5[enb][stb][st1b] + ext_[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({energy, {en + 1, -1, EXT}, {st + 1, en, DP_P}, {st + 1, CTD_5_DANGLE}});

      // .(   ).<   > Terminal mismatch
      energy = base11 + em_->terminal[en1b][enb][stb][st1b] + ext_[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P}, {st + 1, CTD_MISMATCH}});

      if (en < N - 1) {
        // .(   ).<(   ) > Left coax
        energy = base11 + em_->MismatchCoaxial(en1b, enb, stb, st1b);
        if (energy + ext_[en + 1][EXT_GU] <= delta)
          exps.push_back(
              {energy + ext_[en + 1][EXT_GU], {en + 1, -1, EXT_GU}, {st + 1, en - 1, DP_P},
                  {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT}});
        if (energy + ext_[en + 1][EXT_WC] <= delta)
          exps.push_back(
              {energy + ext_[en + 1][EXT_WC], {en + 1, -1, EXT_WC}, {st + 1, en - 1, DP_P},
                  {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT}});
      }

      if (en < N - 2) {
        // (   )<.(   ). > Right coax forward
        energy = base00 + ext_[en + 1][EXT_RC];
        if (energy <= delta)
          exps.push_back({energy, {en + 1, -1, EXT_RC}, {st, en, DP_P}, {en + 2, CTD_RC_WITH_PREV},
              {st, CTD_RC_WITH_NEXT}});
      }

      // (   )(<   ) > Flush coax
      energy = base01 + em_->stack[en1b][enb][WcPair(enb)][stb] + ext_[en][EXT_WC];
      if (energy <= delta)
        exps.push_back({energy, {en, -1, EXT_WC}, {st, en - 1, DP_P}, {en, CTD_FCOAX_WITH_PREV},
            {st, CTD_FCOAX_WITH_NEXT}});

      if (IsGu(enb)) {
        energy = base01 + em_->stack[en1b][enb][GuPair(enb)][stb] + ext_[en][EXT_GU];
        if (energy <= delta)
          exps.push_back({energy, {en, -1, EXT_GU}, {st, en - 1, DP_P}, {en, CTD_FCOAX_WITH_PREV},
              {st, CTD_FCOAX_WITH_NEXT}});
      }
    }
    // Finished exterior loop, don't do anymore.
    return exps;
  }

  // Declare the usual base aliases.
  const auto stb = r_[st];
  const auto st1b = r_[st + 1];
  const auto st2b = r_[st + 2];
  const auto enb = r_[en];
  const auto en1b = r_[en - 1];
  const auto en2b = r_[en - 2];

  // Normal stuff
  if (a == DP_P) {
    // Two loops.
    const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
    for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
      for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
        energy = pc_.TwoLoop(st, en, ist, ien) + dp_[ist][ien][DP_P] - dp_[st][en][a];
        if (energy <= delta) exps.push_back({energy, {ist, ien, DP_P}});
      }
    }

    // Hairpin loop
    energy = pc_.Hairpin(st, en) - dp_[st][en][a];
    if (energy <= delta) exps.emplace_back(energy);

    auto base_and_branch = pc_.augubranch[stb][enb] + em_->multiloop_hack_a - dp_[st][en][a];
    // (<   ><    >)
    energy = base_and_branch + dp_[st + 1][en - 1][DP_U2];
    if (energy <= delta) exps.push_back({energy, {st + 1, en - 1, DP_U2}, {en, CTD_UNUSED}});
    // (3<   ><   >) 3'
    energy = base_and_branch + dp_[st + 2][en - 1][DP_U2] + em_->dangle3[stb][st1b][enb];
    if (energy <= delta) exps.push_back({energy, {st + 2, en - 1, DP_U2}, {en, CTD_3_DANGLE}});
    // (<   ><   >5) 5'
    energy = base_and_branch + dp_[st + 1][en - 2][DP_U2] + em_->dangle5[stb][en1b][enb];
    if (energy <= delta) exps.push_back({energy, {st + 1, en - 2, DP_U2}, {en, CTD_5_DANGLE}});
    // (.<   ><   >.) Terminal mismatch
    energy = base_and_branch + dp_[st + 2][en - 2][DP_U2] + em_->terminal[stb][st1b][en1b][enb];
    if (energy <= delta) exps.push_back({energy, {st + 2, en - 2, DP_U2}, {en, CTD_MISMATCH}});

    for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
      const Base pl1b = r_[piv - 1];
      const Base plb = r_[piv];
      const Base prb = r_[piv + 1];
      const Base pr1b = r_[piv + 2];

      // (.(   )   .) Left outer coax - P
      auto outer_coax = em_->MismatchCoaxial(stb, st1b, en1b, enb);
      energy = base_and_branch + dp_[st + 2][piv][DP_P] + pc_.augubranch[st2b][plb] +
          dp_[piv + 1][en - 2][DP_U] + outer_coax;
      if (energy <= delta)
        exps.push_back({energy, {st + 2, piv, DP_P}, {piv + 1, en - 2, DP_U},
            {st + 2, CTD_LCOAX_WITH_PREV}, {en, CTD_LCOAX_WITH_NEXT}});

      // (.   (   ).) Right outer coax
      energy = base_and_branch + dp_[st + 2][piv][DP_U] + pc_.augubranch[prb][en2b] +
          dp_[piv + 1][en - 2][DP_P] + outer_coax;
      if (energy <= delta)
        exps.push_back({energy, {st + 2, piv, DP_U}, {piv + 1, en - 2, DP_P},
            {piv + 1, CTD_RC_WITH_NEXT}, {en, CTD_RC_WITH_PREV}});

      // (.(   ).   ) Left inner coax
      energy = base_and_branch + dp_[st + 2][piv - 1][DP_P] + pc_.augubranch[st2b][pl1b] +
          dp_[piv + 1][en - 1][DP_U] + em_->MismatchCoaxial(pl1b, plb, st1b, st2b);
      if (energy <= delta)
        exps.push_back({energy, {st + 2, piv - 1, DP_P}, {piv + 1, en - 1, DP_U},
            {st + 2, CTD_RC_WITH_PREV}, {en, CTD_RC_WITH_NEXT}});

      // (   .(   ).) Right inner coax
      energy = base_and_branch + dp_[st + 1][piv][DP_U] + pc_.augubranch[pr1b][en2b] +
          dp_[piv + 2][en - 2][DP_P] + em_->MismatchCoaxial(en2b, en1b, prb, pr1b);
      if (energy <= delta)
        exps.push_back({energy, {st + 1, piv, DP_U}, {piv + 2, en - 2, DP_P},
            {piv + 2, CTD_LCOAX_WITH_NEXT}, {en, CTD_LCOAX_WITH_PREV}});

      // ((   )   ) Left flush coax
      energy = base_and_branch + dp_[st + 1][piv][DP_P] + pc_.augubranch[st1b][plb] +
          dp_[piv + 1][en - 1][DP_U] + em_->stack[stb][st1b][plb][enb];
      if (energy <= delta)
        exps.push_back({energy, {st + 1, piv, DP_P}, {piv + 1, en - 1, DP_U},
            {st + 1, CTD_FCOAX_WITH_PREV}, {en, CTD_FCOAX_WITH_NEXT}});

      // (   (   )) Right flush coax
      energy = base_and_branch + dp_[st + 1][piv][DP_U] + pc_.augubranch[prb][en1b] +
          dp_[piv + 1][en - 1][DP_P] + em_->stack[stb][prb][en1b][enb];
      if (energy <= delta)
        exps.push_back({energy, {st + 1, piv, DP_U}, {piv + 1, en - 1, DP_P},
            {piv + 1, CTD_FCOAX_WITH_NEXT}, {en, CTD_FCOAX_WITH_PREV}});
    }
    return exps;
  }

  // Left unpaired. Either DP_U or DP_U2.
  if (st + 1 < en && (a == DP_U || a == DP_U2)) {
    energy = dp_[st + 1][en][a] - dp_[st][en][a];
    if (energy <= delta) exps.push_back({energy, {st + 1, en, a}});
  }

  // Pair here.
  for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
    //   (   .   )<   (
    // stb pl1b pb   pr1b
    auto pb = r_[piv];
    auto pl1b = r_[piv - 1];
    // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
    auto base00 = dp_[st][piv][DP_P] + pc_.augubranch[stb][pb] - dp_[st][en][a];
    auto base01 = dp_[st][piv - 1][DP_P] + pc_.augubranch[stb][pl1b] - dp_[st][en][a];
    auto base10 = dp_[st + 1][piv][DP_P] + pc_.augubranch[st1b][pb] - dp_[st][en][a];
    auto base11 = dp_[st + 1][piv - 1][DP_P] + pc_.augubranch[st1b][pl1b] - dp_[st][en][a];

    // Check a == U_RC:
    // (   )<.( ** ). > Right coax backward
    if (a == DP_U_RC) {
      energy = base11 + em_->MismatchCoaxial(pl1b, pb, stb, st1b);
      // Our ctds will have already been set by now.
      if (energy <= delta) exps.push_back({energy, {st + 1, piv - 1, DP_P}});
      if (energy + dp_[piv + 1][en][DP_U] <= delta)
        exps.push_back(
            {energy + dp_[piv + 1][en][DP_U], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U}});
      continue;
    }
    // From here on, a must be U, U2, U_WC, or U_GU.

    // (   )<   > - U, U2, U_WC?, U_GU?
    energy = base00;
    if (a == DP_U) {
      if (energy <= delta) exps.push_back({energy, {st, piv, DP_P}, {st, CTD_UNUSED}});
      if (energy + dp_[piv + 1][en][DP_U] <= delta)
        exps.push_back({energy + dp_[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U},
            {st, CTD_UNUSED}});
    }
    if (a == DP_U2 && energy + dp_[piv + 1][en][DP_U] <= delta)
      exps.push_back({energy + dp_[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U},
          {st, CTD_UNUSED}});
    if (a == DP_U_WC || a == DP_U_GU) {
      // Make sure we don't form any branches that are not the right type of pair.
      if ((a == DP_U_WC && IsWcPair(stb, pb)) || (a == DP_U_GU && IsGuPair(stb, pb))) {
        if (energy <= delta) exps.push_back({energy, {st, piv, DP_P}});
        if (energy + dp_[piv + 1][en][DP_U] <= delta)
          exps.push_back({energy + dp_[piv + 1][en][DP_U], {st, piv, DP_P}, {piv + 1, en, DP_U}});
      }
      continue;
    }
    // From here on, a must be U or U2.
    assert(a == DP_U || a == DP_U2);

    // (   )3<   > 3' - U, U2
    energy = base01 + em_->dangle3[pl1b][pb][stb];
    // Can only let the rest be unpaired if we only need one branch, i.e. DP_U not DP_U2.
    if (a == DP_U && energy <= delta)
      exps.push_back({energy, {st, piv - 1, DP_P}, {st, CTD_3_DANGLE}});
    if (energy + dp_[piv + 1][en][DP_U] <= delta)
      exps.push_back({energy + dp_[piv + 1][en][DP_U], {st, piv - 1, DP_P}, {piv + 1, en, DP_U},
          {st, CTD_3_DANGLE}});

    // 5(   )<   > 5' - U, U2
    energy = base10 + em_->dangle5[pb][stb][st1b];
    if (a == DP_U && energy <= delta)
      exps.push_back({energy, {st + 1, piv, DP_P}, {st + 1, CTD_5_DANGLE}});
    if (energy + dp_[piv + 1][en][DP_U] <= delta)
      exps.push_back({energy + dp_[piv + 1][en][DP_U], {st + 1, piv, DP_P}, {piv + 1, en, DP_U},
          {st + 1, CTD_5_DANGLE}});

    // .(   ).<   > Terminal mismatch - U, U2
    energy = base11 + em_->terminal[pl1b][pb][stb][st1b];
    if (a == DP_U && energy <= delta)
      exps.push_back({energy, {st + 1, piv - 1, DP_P}, {}, {st + 1, CTD_MISMATCH}});
    if (energy + dp_[piv + 1][en][DP_U] <= delta)
      exps.push_back({energy + dp_[piv + 1][en][DP_U], {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U},
          {st + 1, CTD_MISMATCH}});

    // .(   ).<(   ) > Left coax - U, U2
    energy = base11 + em_->MismatchCoaxial(pl1b, pb, stb, st1b);
    if (energy + dp_[piv + 1][en][DP_U_WC] <= delta)
      exps.push_back({energy + dp_[piv + 1][en][DP_U_WC], {st + 1, piv - 1, DP_P},
          {piv + 1, en, DP_U_WC}, {st + 1, CTD_LCOAX_WITH_NEXT}, {piv + 1, CTD_LCOAX_WITH_PREV}});
    if (energy + dp_[piv + 1][en][DP_U_GU] <= delta)
      exps.push_back({energy + dp_[piv + 1][en][DP_U_GU], {st + 1, piv - 1, DP_P},
          {piv + 1, en, DP_U_GU}, {st + 1, CTD_LCOAX_WITH_NEXT}, {piv + 1, CTD_LCOAX_WITH_PREV}});

    // (   )<.(   ). > Right coax forward - U, U2
    energy = base00 + dp_[piv + 1][en][DP_U_RC];
    if (energy <= delta)
      exps.push_back({energy, {st, piv, DP_P}, {piv + 1, en, DP_U_RC}, {st, CTD_RC_WITH_NEXT},
          {piv + 2, CTD_RC_WITH_PREV}});

    // (   )(<   ) > Flush coax - U, U2
    energy = base01 + em_->stack[pl1b][pb][WcPair(pb)][stb] + dp_[piv][en][DP_U_WC];
    if (energy <= delta)
      exps.push_back({energy, {st, piv - 1, DP_P}, {piv, en, DP_U_WC}, {st, CTD_FCOAX_WITH_NEXT},
          {piv, CTD_FCOAX_WITH_PREV}});

    if (IsGu(pb)) {
      energy = base01 + em_->stack[pl1b][pb][GuPair(pb)][stb] + dp_[piv][en][DP_U_GU];
      if (energy <= delta)
        exps.push_back({energy, {st, piv - 1, DP_P}, {piv, en, DP_U_GU}, {st, CTD_FCOAX_WITH_NEXT},
            {piv, CTD_FCOAX_WITH_PREV}});
    }
  }

  return exps;
}

}  // namespace mrna::subopt::t04
