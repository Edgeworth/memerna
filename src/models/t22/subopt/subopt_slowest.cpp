// Copyright 2023 Eliot Courtney.
#include "models/t22/subopt/subopt_slowest.h"

#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <exception>
#include <memory>
#include <utility>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "api/subopt/subopt.h"
#include "api/trace/trace.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "util/error.h"

namespace mrna::md::t22 {

using t04::DP_P;
using t04::EXT;

SuboptSlowest::SuboptSlowest(Primary r, Model::Ptr em, DpState dp, SuboptCfg cfg)
    : r_(std::move(r)), em_(std::move(em)), dp_(std::move(dp)), cfg_(cfg) {}

int SuboptSlowest::Run(const SuboptCallback& fn) {
  res_ = SuboptResult(ZERO_E, trace::TraceResult(Secondary(r_.size()), Ctds(r_.size())));
  q_.reserve(r_.size());
  cache_.Reserve(r_.size());

  verify(em_->cfg.lonely_pairs != erg::EnergyCfg::LonelyPairs::OFF,
      "fully disallowing lonely pairs is not supported in this energy model");
  verify(em_->cfg.ctd == erg::EnergyCfg::Ctd::ALL,
      "only full CTDs are supported in this energy model");

  spdlog::debug("t22 {} with cfg {}", __func__, em_->cfg);

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

std::pair<int, Energy> SuboptSlowest::RunInternal(
    const SuboptCallback& fn, Energy delta, bool exact_energy, int max) {
  int count = 0;
  Energy next_seen = MAX_E;
  Energy energy = ZERO_E;
  Energy mfe = dp_.t04.ext[0][EXT];
  q_.clear();
  unexpanded_.clear();
  q_.push_back({.child_idx = 0, .to_expand = t04::DpIndex(0, -1, EXT), .should_unexpand = false});
  while (!q_.empty()) {
    auto& s = q_.back();
    assert(s.to_expand.has_value());
    auto to_expand = *s.to_expand;

    const auto& exps = GetExpansion(to_expand);
    assert(!exps.empty());

    // Go to next child:
    if (s.child_idx != 0) {
      const auto& pexp = exps[s.child_idx - 1];
      pexp.ctd0.MaybeRemove(res_.tb.ctd);
      pexp.ctd1.MaybeRemove(res_.tb.ctd);
      pexp.pair.MaybeRemove(res_.tb.s);
      if (pexp.idx1.has_value()) unexpanded_.pop_back();
      energy -= pexp.delta;
    }

    if (s.child_idx != static_cast<int>(exps.size()) && exps[s.child_idx].delta + energy > delta)
      next_seen = std::min(next_seen, exps[s.child_idx].delta + energy);

    if (s.child_idx == static_cast<int>(exps.size()) || exps[s.child_idx].delta + energy > delta) {
      if (s.should_unexpand) unexpanded_.push_back(to_expand);
      q_.pop_back();
      continue;
    }

    const auto& exp = exps[s.child_idx++];
    DfsState ns = {0, exp.idx0, false};
    energy += exp.delta;

    if (!exp.idx0.has_value()) {
      assert(!exp.idx1.has_value());
      assert(!exp.ctd0.IsValid() && !exp.ctd1.IsValid());
      assert(!exp.pair.IsValid());

      if (unexpanded_.empty()) {
        // At a terminal state.
        if (!exact_energy || energy == delta) {
          res_.energy = energy + mfe;
          fn(res_);
          ++count;

          if (count == max) return {count, CAP_E};
        }
        continue;
      }

      ns.to_expand = unexpanded_.back();
      unexpanded_.pop_back();
      ns.should_unexpand = true;
    } else {
      // Update global state with this expansion.
      exp.ctd0.MaybeApply(res_.tb.ctd);
      exp.ctd1.MaybeApply(res_.tb.ctd);
      exp.pair.MaybeApply(res_.tb.s);
      if (exp.idx1.has_value()) unexpanded_.push_back(*exp.idx1);
    }

    q_.push_back(ns);
  }
  assert(unexpanded_.empty() && energy == ZERO_E && res_.tb.s == Secondary(res_.tb.s.size()) &&
      res_.tb.ctd == Ctds(res_.tb.ctd.size()));
  return {count, next_seen};
}

std::vector<Expansion> SuboptSlowest::GenerateExpansions(
    const DpIndex& to_expand, Energy delta) const {
  if (std::holds_alternative<t04::DpIndex>(to_expand)) {
    auto idx = std::get<t04::DpIndex>(to_expand);
    int st = idx.st;
    int en = idx.en;
    int a = idx.a;
    if (en == -1) {
      return ExtExpansions(st, en, a, delta);
    } else {
      if (a == DP_P) {
        return PairedOrNoStackExpansions(st, en, /*is_nostack=*/false, delta);
      } else {
        return UnpairedExpansions(st, en, a, delta);
      }
    }
  } else if (std::holds_alternative<PenultimateIndex>(to_expand)) {
    auto idx = std::get<PenultimateIndex>(to_expand);
    const int st = idx.st;
    const int en = idx.en;
    const int length = idx.len;
    return PenultimateExpansions(st, en, length, delta);
  } else {
    auto idx = std::get<NoStackIndex>(to_expand);
    const int st = idx.st;
    const int en = idx.en;
    return PairedOrNoStackExpansions(st, en, /*is_nostack=*/true, delta);
  }
  // Should not reach here.
  bug();
  return {};
}

std::vector<Expansion> SuboptSlowest::ExtExpansions(int st, int en, int a, Energy delta) const {
  // Case: No pair starting here
  if (a == EXT && st + 1 < N && ext[st + 1][EXT] + em.PfUnpaired(st) == ext[st][EXT]) {
    next.push_back({.idx0 = t04::DpIndex(st + 1, -1, EXT)});
  }
  for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
    // .   .   .   (   .   .   .   )   <   >
    //           stb  st1b   en1b  enb   rem
    const auto stb = r[st];
    const auto st1b = r[st + 1];
    const auto enb = r[en];
    const auto en1b = r[en - 1];
    const auto base00 = dp[st][en][DP_P] + em.AuGuPenalty(stb, enb);
    const auto base01 = dp[st][en - 1][DP_P] + em.AuGuPenalty(stb, en1b);
    const auto base10 = dp[st + 1][en][DP_P] + em.AuGuPenalty(st1b, enb);
    const auto base11 = dp[st + 1][en - 1][DP_P] + em.AuGuPenalty(st1b, en1b);

    // (   )<.( * ). > Right coax backward
    if (a == EXT_RC && em.cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
      // Don't set CTDs here since they will have already been set.
      if (base11 + em.MismatchCoaxial(en1b, enb, stb, st1b) + em.PfUnpaired(st) +
              em.PfUnpaired(en) + ext[en + 1][EXT] ==
          ext[st][EXT_RC]) {
        next.push_back(
            {.idx0 = t04::DpIndex(st + 1, en - 1, DP_P), .idx1 = t04::DpIndex(en + 1, -1, EXT)});
      }
    }

    // (   )<   >
    auto val = base00 + ext[en + 1][EXT];
    if (val == ext[st][a] && (a != EXT_WC || IsWcPair(stb, enb)) &&
        (a != EXT_GU || IsGuPair(stb, enb))) {
      // EXT_WC and EXT_GU will have already had their ctds set.
      Expansion exp{.idx0 = t04::DpIndex(st, en, DP_P), .idx1 = t04::DpIndex(en + 1, -1, EXT)};
      if (a == EXT) exp.ctd0 = {st, CTD_UNUSED};
      next.push_back(exp);
    }

    // Only look at EXT from here on.
    if (a != EXT) continue;

    // (   )3<   > 3'
    if (base01 + em.dangle3[en1b][enb][stb] + em.PfUnpaired(en) + ext[en + 1][EXT] ==
        ext[st][EXT]) {
      next.push_back({
          .idx0 = t04::DpIndex(st, en - 1, DP_P),
          .idx1 = t04::DpIndex(en + 1, -1, EXT),
          .ctd0{st, CTD_3_DANGLE},
      });
    }
    // 5(   )<   > 5'
    if (base10 + em.dangle5[enb][stb][st1b] + em.PfUnpaired(st) + ext[en + 1][EXT] ==
        ext[st][EXT]) {
      next.push_back({.idx0 = t04::DpIndex(st + 1, en, DP_P),
          .idx1 = t04::DpIndex(en + 1, -1, EXT),
          .ctd0{st + 1, CTD_5_DANGLE}});
    }
    // .(   ).<   > Terminal mismatch
    if (base11 + em.terminal[en1b][enb][stb][st1b] + em.PfUnpaired(st) + em.PfUnpaired(en) +
            ext[en + 1][EXT] ==
        ext[st][EXT]) {
      next.push_back({.idx0 = t04::DpIndex(st + 1, en - 1, DP_P),
          .idx1 = t04::DpIndex(en + 1, -1, EXT),
          .ctd0{st + 1, CTD_MISMATCH}});
    }

    if (en < N - 1 && em.cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
      // .(   ).<(   ) > Left coax  x
      val =
          base11 + em.MismatchCoaxial(en1b, enb, stb, st1b) + em.PfUnpaired(st) + em.PfUnpaired(en);
      if (val + ext[en + 1][EXT_WC] == ext[st][EXT]) {
        next.push_back({.idx0 = t04::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT_WC),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{en + 1, CTD_LCOAX_WITH_PREV}});
      }
      if (val + ext[en + 1][EXT_GU] == ext[st][EXT]) {
        next.push_back({.idx0 = t04::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT_GU),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{en + 1, CTD_LCOAX_WITH_PREV}});
      }

      // (   )<.(   ). > Right coax forward
      if (base00 + ext[en + 1][EXT_RC] == ext[st][EXT]) {
        next.push_back({.idx0 = t04::DpIndex(st, en, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT_RC),
            .ctd0{st, CTD_RC_WITH_NEXT},
            .ctd1{en + 2, CTD_RC_WITH_PREV}});
      }

      // (   )(<   ) > Flush coax
      if (base01 + em.stack[en1b][enb][WcPair(enb)][stb] + ext[en][EXT_WC] == ext[st][EXT]) {
        next.push_back({.idx0 = t04::DpIndex(st, en - 1, DP_P),
            .idx1 = t04::DpIndex(en, -1, EXT_WC),
            .ctd0{st, CTD_FCOAX_WITH_NEXT},
            .ctd1{en, CTD_FCOAX_WITH_PREV}});
      }
      if (IsGu(enb) &&
          base01 + em.stack[en1b][enb][GuPair(enb)][stb] + ext[en][EXT_GU] == ext[st][EXT]) {
        next.push_back({.idx0 = t04::DpIndex(st, en - 1, DP_P),
            .idx1 = t04::DpIndex(en, -1, EXT_GU),
            .ctd0{st, CTD_FCOAX_WITH_NEXT},
            .ctd1{en, CTD_FCOAX_WITH_PREV}});
      }
    }
  }
}

std::vector<Expansion> SuboptSlowest::PairedOrNoStackExpansions(
    int st, int en, bool is_nostack, Energy delta) const {
  const auto stb = r[st];
  const auto st1b = r[st + 1];
  const auto st2b = r[st + 2];
  const auto enb = r[en];
  const auto en1b = r[en - 1];
  const auto en2b = r[en - 2];

  if (!is_nostack) {
    const int max_stack = en - st - HAIRPIN_MIN_SZ + 1;
    const Energy bulge_left = em.Bulge(r, st, en, st + 2, en - 1);
    const Energy bulge_right = em.Bulge(r, st, en, st + 1, en - 2);

    for (int length = 2; 2 * length <= max_stack; ++length) {
      auto none = em.stack[r[st]][r[st + 1]][r[en - 1]][r[en]] +
          em.penultimate_stack[en1b][enb][stb][st1b] + em.PfPaired(st, en);
      if (length == 2 &&
          none + nostack[st + 1][en - 1] +
                  em.penultimate_stack[r[st]][r[st + 1]][r[en - 1]][r[en]] ==
              dp[st][en][DP_P]) {
        next.push_back({.idx0 = NoStackIndex(st + 1, en - 1), .pair{st, en}});
      }
      if (none + penult[st + 1][en - 1][length - 1] == dp[st][en][DP_P]) {
        next.push_back({.idx0 = PenultimateIndex(st + 1, en - 1, length - 1), .pair{st, en}});
      }

      auto left = bulge_left + em.penultimate_stack[en1b][enb][stb][st2b];
      if (length == 2 &&
          left + nostack[st + 2][en - 1] +
                  em.penultimate_stack[r[st]][r[st + 2]][r[en - 1]][r[en]] ==
              dp[st][en][DP_P]) {
        next.push_back({.idx0 = NoStackIndex(st + 2, en - 1), .pair{st, en}});
      }
      if (left + penult[st + 2][en - 1][length - 1] == dp[st][en][DP_P]) {
        next.push_back({.idx0 = PenultimateIndex(st + 2, en - 1, length - 1), .pair{st, en}});
      }

      auto right = bulge_right + em.penultimate_stack[en2b][enb][stb][st1b];
      if (length == 2 &&
          right + nostack[st + 1][en - 2] +
                  em.penultimate_stack[r[st]][r[st + 1]][r[en - 2]][r[en]] ==
              dp[st][en][DP_P]) {
        next.push_back({.idx0 = NoStackIndex(st + 1, en - 2), .pair{st, en}});
      }
      if (right + penult[st + 1][en - 2][length - 1] == dp[st][en][DP_P]) {
        next.push_back({.idx0 = PenultimateIndex(st + 1, en - 2, length - 1), .pair{st, en}});
      }
    }
  }

  const auto target = is_nostack ? nostack[st][en] : dp[st][en][DP_P];

  // Following largely matches the above DP so look up there for comments.
  const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
  for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
    for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
      // Try all internal loops. We don't check stacks or 1 nuc bulge loops.
      if (dp[ist][ien][DP_P] < CAP_E && ist - st + en - ien > 3) {
        const auto val = em.TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P];
        if (val == target) next.push_back({.idx0 = t04::DpIndex(ist, ien, DP_P), .pair{st, en}});
      }
    }
  }

  if (em.Hairpin(r, st, en) == target) {
    next.push_back({.pair{st, en}});
  }

  const auto base_branch_cost =
      em.AuGuPenalty(stb, enb) + em.PfPaired(st, en) + em.multiloop_hack_a + em.multiloop_hack_b;
  // (<   ><    >)
  if (base_branch_cost + dp[st + 1][en - 1][DP_U2] == target) {
    next.push_back({
        .idx0 = t04::DpIndex(st + 1, en - 1, DP_U2),
        .ctd0{en, CTD_UNUSED},
        .pair{st, en},
    });
  }
  // (3<   ><   >) 3'
  if (base_branch_cost + dp[st + 2][en - 1][DP_U2] + em.dangle3[stb][st1b][enb] +
          em.PfUnpaired(st + 1) ==
      target) {
    next.push_back(
        {.idx0 = t04::DpIndex(st + 2, en - 1, DP_U2), .ctd0{en, CTD_3_DANGLE}, .pair{st, en}});
  }
  // (<   ><   >5) 5'
  if (base_branch_cost + dp[st + 1][en - 2][DP_U2] + em.dangle5[stb][en1b][enb] +
          em.PfUnpaired(en - 1) ==
      target) {
    next.push_back(
        {.idx0 = t04::DpIndex(st + 1, en - 2, DP_U2), .ctd0{en, CTD_5_DANGLE}, .pair{st, en}});
  }
  // (.<   ><   >.) Terminal mismatch
  if (base_branch_cost + dp[st + 2][en - 2][DP_U2] + em.terminal[stb][st1b][en1b][enb] +
          em.PfUnpaired(st + 1) + em.PfUnpaired(en - 1) ==
      target) {
    next.push_back(
        {.idx0 = t04::DpIndex(st + 2, en - 2, DP_U2), .ctd0{en, CTD_MISMATCH}, .pair{st, en}});
  }

  if (em.cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
    for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
      const Base pl1b = r[piv - 1];
      const Base plb = r[piv];
      const Base prb = r[piv + 1];
      const Base pr1b = r[piv + 2];

      // (.(   )   .) Left outer coax - P
      const auto outer_coax =
          em.MismatchCoaxial(stb, st1b, en1b, enb) + em.PfUnpaired(st + 1) + em.PfUnpaired(en - 1);
      if (base_branch_cost + dp[st + 2][piv][DP_P] + em.multiloop_hack_b +
              em.AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax ==
          target) {
        next.push_back({.idx0 = t04::DpIndex(st + 2, piv, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en - 2, DP_U),
            .ctd0{en, CTD_LCOAX_WITH_NEXT},
            .ctd1{st + 2, CTD_LCOAX_WITH_PREV},
            .pair{st, en}});
      }
      // (.   (   ).) Right outer coax
      if (base_branch_cost + dp[st + 2][piv][DP_U] + em.multiloop_hack_b +
              em.AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax ==
          target) {
        next.push_back({.idx0 = t04::DpIndex(st + 2, piv, DP_U),
            .idx1 = t04::DpIndex(piv + 1, en - 2, DP_P),
            .ctd0{en, CTD_RC_WITH_PREV},
            .ctd1{piv + 1, CTD_RC_WITH_NEXT},
            .pair{st, en}});
      }

      // (.(   ).   ) Left inner coax
      if (base_branch_cost + dp[st + 2][piv - 1][DP_P] + em.multiloop_hack_b +
              em.AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
              em.MismatchCoaxial(pl1b, plb, st1b, st2b) + em.PfUnpaired(st + 1) +
              em.PfUnpaired(piv) ==
          target) {
        next.push_back({.idx0 = t04::DpIndex(st + 2, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en - 1, DP_U),
            .ctd0{en, CTD_RC_WITH_NEXT},
            .ctd1{st + 2, CTD_RC_WITH_PREV},
            .pair{st, en}});
      }
      // (   .(   ).) Right inner coax
      if (base_branch_cost + dp[st + 1][piv][DP_U] + em.multiloop_hack_b +
              em.AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
              em.MismatchCoaxial(en2b, en1b, prb, pr1b) + em.PfUnpaired(piv + 1) +
              em.PfUnpaired(en - 1) ==
          target) {
        next.push_back({.idx0 = t04::DpIndex(st + 1, piv, DP_U),
            .idx1 = t04::DpIndex(piv + 2, en - 2, DP_P),
            .ctd0{en, CTD_LCOAX_WITH_PREV},
            .ctd1{piv + 2, CTD_LCOAX_WITH_NEXT},
            .pair{st, en}});
      }

      // ((   )   ) Left flush coax
      if (base_branch_cost + dp[st + 1][piv][DP_P] + em.multiloop_hack_b +
              em.AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] +
              em.stack[stb][st1b][plb][enb] ==
          target) {
        next.push_back({.idx0 = t04::DpIndex(st + 1, piv, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en - 1, DP_U),
            .ctd0{en, CTD_FCOAX_WITH_NEXT},
            .ctd1{st + 1, CTD_FCOAX_WITH_PREV},
            .pair{st, en}});
      }
      // (   (   )) Right flush coax
      if (base_branch_cost + dp[st + 1][piv][DP_U] + em.multiloop_hack_b +
              em.AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] +
              em.stack[stb][prb][en1b][enb] ==
          target) {
        next.push_back({.idx0 = t04::DpIndex(st + 1, piv, DP_U),
            .idx1 = t04::DpIndex(piv + 1, en - 1, DP_P),
            .ctd0{en, CTD_FCOAX_WITH_PREV},
            .ctd1{piv + 1, CTD_FCOAX_WITH_NEXT},
            .pair{st, en}});
      }
    }
  }
}

std::vector<Expansion> SuboptSlowest::UnpairedExpansions(
    int st, int en, int a, Energy delta) const {
  const auto stb = r[st];
  const auto st1b = r[st + 1];

  // Deal with the rest of the cases:
  // Left unpaired. Either DP_U or DP_U2.
  if (st + 1 < en && (a == DP_U || a == DP_U2) &&
      dp[st + 1][en][a] + em.PfUnpaired(st) == dp[st][en][a]) {
    next.push_back({.idx0 = t04::DpIndex(st + 1, en, a)});
  }

  // Pair here.
  for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
    //   (   .   )<   (
    // stb pl1b pb   pr1b
    const auto pb = r[piv];
    const auto pl1b = r[piv - 1];
    // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the
    // right.
    const auto base00 = dp[st][piv][DP_P] + em.AuGuPenalty(stb, pb) + em.multiloop_hack_b;
    const auto base01 = dp[st][piv - 1][DP_P] + em.AuGuPenalty(stb, pl1b) + em.multiloop_hack_b;
    const auto base10 = dp[st + 1][piv][DP_P] + em.AuGuPenalty(st1b, pb) + em.multiloop_hack_b;
    const auto base11 =
        dp[st + 1][piv - 1][DP_P] + em.AuGuPenalty(st1b, pl1b) + em.multiloop_hack_b;

    // Min is for either placing another unpaired or leaving it as nothing.
    // If we're at U2, don't allow leaving as nothing.
    auto right_unpaired = dp[piv + 1][en][DP_U];

    // We need to store both of these because in the case that, due to
    // pseudofree energies, both are equal, we need to consider both
    // possibilities for a randomised traceback.
    bool can_right_paired = true;
    bool can_right_unpaired = false;
    if (a != DP_U2) {
      const auto unpaired_cum = em.PfUnpairedCum(piv + 1, en);
      if (unpaired_cum == right_unpaired) can_right_unpaired = true;
      if (unpaired_cum < right_unpaired) {
        can_right_paired = false;
        can_right_unpaired = true;
        right_unpaired = unpaired_cum;
      }
    }

    if (em.cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
      // Check a == U_RC:
      // (   )<.( ** ). > Right coax backward
      if (a == DP_U_RC) {
        if (base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b) + em.PfUnpaired(st) +
                em.PfUnpaired(piv) + right_unpaired ==
            dp[st][en][DP_U_RC]) {
          // Ctds were already set from the recurrence that called this.
          Expansion exp{.idx0 = t04::DpIndex(st + 1, piv - 1, DP_P)};
          if (can_right_unpaired) next.push_back(exp);
          if (can_right_paired) {
            exp.idx1 = t04::DpIndex(piv + 1, en, DP_U);
            next.push_back(exp);
          }
        }
      }
    }

    // (   )<   > - U, U2, U_WC?, U_GU?
    if (base00 + right_unpaired == dp[st][en][a] && (a != DP_U_WC || IsWcPair(stb, pb)) &&
        (a != DP_U_GU || IsGuPair(stb, pb))) {
      // If U_WC, or U_GU, we were involved in some sort of coaxial stack previously, and
      // were already set.
      Expansion exp{.idx0 = t04::DpIndex(st, piv, DP_P)};
      if (a != DP_U_WC && a != DP_U_GU) exp.ctd0 = {st, CTD_UNUSED};
      if (can_right_unpaired) next.push_back(exp);
      if (can_right_paired) {
        exp.idx1 = t04::DpIndex(piv + 1, en, DP_U);
        next.push_back(exp);
      }
    }

    // The rest of the cases are for U and U2.
    if (a != DP_U && a != DP_U2) continue;

    // (   )3<   > 3' - U, U2
    if (base01 + em.dangle3[pl1b][pb][stb] + em.PfUnpaired(piv) + right_unpaired == dp[st][en][a]) {
      Expansion exp{.idx0 = t04::DpIndex(st, piv - 1, DP_P), .ctd0{st, CTD_3_DANGLE}};
      if (can_right_unpaired) next.push_back(exp);
      if (can_right_paired) {
        exp.idx1 = t04::DpIndex(piv + 1, en, DP_U);
        next.push_back(exp);
      }
    }
    // 5(   )<   > 5' - U, U2
    if (base10 + em.dangle5[pb][stb][st1b] + em.PfUnpaired(st) + right_unpaired == dp[st][en][a]) {
      Expansion exp{.idx0 = t04::DpIndex(st + 1, piv, DP_P), .ctd0{st + 1, CTD_5_DANGLE}};
      if (can_right_unpaired) next.push_back(exp);
      if (can_right_paired) {
        exp.idx1 = t04::DpIndex(piv + 1, en, DP_U);
        next.push_back(exp);
      }
    }
    // .(   ).<   > Terminal mismatch - U, U2
    if (base11 + em.terminal[pl1b][pb][stb][st1b] + em.PfUnpaired(st) + em.PfUnpaired(piv) +
            right_unpaired ==
        dp[st][en][a]) {
      Expansion exp{.idx0 = t04::DpIndex(st + 1, piv - 1, DP_P), .ctd0{st + 1, CTD_MISMATCH}};
      if (can_right_unpaired) next.push_back(exp);
      if (can_right_paired) {
        exp.idx1 = t04::DpIndex(piv + 1, en, DP_U);
        next.push_back(exp);
      }
    }

    if (em.cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
      // .(   ).<(   ) > Left coax - U, U2
      auto val =
          base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b) + em.PfUnpaired(st) + em.PfUnpaired(piv);
      if (val + dp[piv + 1][en][DP_U_WC] == dp[st][en][a]) {
        next.push_back({
            .idx0 = t04::DpIndex(st + 1, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U_WC),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{piv + 1, CTD_LCOAX_WITH_PREV},
        });
      }
      if (val + dp[piv + 1][en][DP_U_GU] == dp[st][en][a]) {
        next.push_back({
            .idx0 = t04::DpIndex(st + 1, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U_GU),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{piv + 1, CTD_LCOAX_WITH_PREV},
        });
      }

      // (   )<.(   ). > Right coax forward - U, U2
      if (base00 + dp[piv + 1][en][DP_U_RC] == dp[st][en][a]) {
        next.push_back({
            .idx0 = t04::DpIndex(st, piv, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U_RC),
            .ctd0{st, CTD_RC_WITH_NEXT},
            .ctd1{piv + 2, CTD_RC_WITH_PREV},
        });
      }

      // (   )(<   ) > Flush coax - U, U2
      if (base01 + em.stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC] == dp[st][en][a]) {
        next.push_back({
            .idx0 = t04::DpIndex(st, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv, en, DP_U_WC),
            .ctd0{st, CTD_FCOAX_WITH_NEXT},
            .ctd1{piv, CTD_FCOAX_WITH_PREV},
        });
      }
      if ((IsGu(pb)) &&
          base01 + em.stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU] == dp[st][en][a]) {
        next.push_back({
            .idx0 = t04::DpIndex(st, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv, en, DP_U_GU),
            .ctd0{st, CTD_FCOAX_WITH_NEXT},
            .ctd1{piv, CTD_FCOAX_WITH_PREV},
        });
      }
    }
  }
}

std::vector<Expansion> SuboptSlowest::PenultimateExpansions(
    int st, int en, int length, Energy delta) const {
  const auto bulge_left = em.Bulge(r, st, en, st + 2, en - 1);
  const auto bulge_right = em.Bulge(r, st, en, st + 1, en - 2);

  auto none = em.stack[r[st]][r[st + 1]][r[en - 1]][r[en]] + em.PfPaired(st, en);
  if (length == 2 &&
      none + nostack[st + 1][en - 1] + em.penultimate_stack[r[st]][r[st + 1]][r[en - 1]][r[en]] ==
          penult[st][en][length]) {
    next.push_back({.idx0 = NoStackIndex(st + 1, en - 1), .pair{st, en}});
  }
  if (none + penult[st + 1][en - 1][length - 1] == penult[st][en][length]) {
    next.push_back({.idx0 = PenultimateIndex(st + 1, en - 1, length - 1), .pair{st, en}});
  }

  auto left = bulge_left;
  if (length == 2 &&
      left + nostack[st + 2][en - 1] + em.penultimate_stack[r[st]][r[st + 2]][r[en - 1]][r[en]] ==
          penult[st][en][length]) {
    next.push_back({.idx0 = NoStackIndex(st + 2, en - 1), .pair{st, en}});
  }
  if (left + penult[st + 2][en - 1][length - 1] == penult[st][en][length]) {
    next.push_back({.idx0 = PenultimateIndex(st + 2, en - 1, length - 1), .pair{st, en}});
  }

  auto right = bulge_right;
  if (length == 2 &&
      right + nostack[st + 1][en - 2] + em.penultimate_stack[r[st]][r[st + 1]][r[en - 2]][r[en]] ==
          penult[st][en][length]) {
    next.push_back({.idx0 = NoStackIndex(st + 1, en - 2), .pair{st, en}});
  }
  if (right + penult[st + 1][en - 2][length - 1] == penult[st][en][length]) {
    next.push_back({.idx0 = PenultimateIndex(st + 1, en - 2, length - 1), .pair{st, en}});
  }
}

}  // namespace mrna::md::t22
