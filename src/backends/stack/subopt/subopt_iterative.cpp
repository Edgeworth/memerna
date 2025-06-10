// Copyright 2023 Eliot Courtney.
#include "backends/stack/subopt/subopt_iterative.h"

#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <utility>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "api/trace/trace.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "util/error.h"

namespace mrna::md::stack {

using base::DP_P;
using base::DP_U;
using base::DP_U2;
using base::DP_U_GU;
using base::DP_U_RC;
using base::DP_U_WC;
using base::EXT;
using base::EXT_GU;
using base::EXT_RC;
using base::EXT_WC;

template <bool UseLru>
SuboptIterative<UseLru>::SuboptIterative(Primary r, Model::Ptr m, DpState dp, SuboptCfg cfg)
    : r_(std::move(r)), m_(std::move(m)), dp_(std::move(dp)), cfg_(cfg),
      cache_(r_, MaxLinearIndex(r_.size())) {}

template <bool UseLru>
int SuboptIterative<UseLru>::Run(const SuboptCallback& fn) {
  res_ = SuboptResult(ZERO_E, trace::TraceResult(Secondary(r_.size()), Ctds(r_.size())));
  q_.reserve(r_.size());

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false, true},
      .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::NONE},
  };
  support.VerifySupported(funcname(), m_->cfg());

  spdlog::debug("stack {} with cfg {}", funcname(), m_->cfg());

  if (cfg_.sorted || cfg_.strucs != SuboptCfg::MAX_STRUCTURES || cfg_.time_secs >= 0.0) {
    int count = 0;
    Energy delta = ZERO_E;
    auto start_time = std::chrono::steady_clock::now();
    while (count < cfg_.strucs && delta != MAX_E && delta <= cfg_.delta) {
      if (cfg_.time_secs >= 0.0) {
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
            std::chrono::steady_clock::now() - start_time);
        if (elapsed.count() >= cfg_.time_secs) break;
      }
      auto res = RunInternal(fn, delta, true, cfg_.strucs - count);
      count += res.first;
      delta = res.second;
    }
    return count;
  }
  return RunInternal(fn, cfg_.delta, false, cfg_.strucs).first;
}

template <bool UseLru>
std::pair<int, Energy> SuboptIterative<UseLru>::RunInternal(
    const SuboptCallback& fn, Energy delta, bool exact_energy, int max) {
  int count = 0;
  Energy next_seen = MAX_E;
  Energy energy = ZERO_E;
  Energy mfe = dp_.base.ext[0][EXT];
  q_.clear();
  unexpanded_.clear();
  q_.push_back({.expand_idx = 0, .to_expand = base::DpIndex(0, -1, EXT), .should_unexpand = false});
  while (!q_.empty()) {
    auto& s = q_.back();
    assert(s.to_expand.has_value());
    auto to_expand = *s.to_expand;

    const auto& exps = GetExpansion(to_expand);
    assert(!exps.empty());

    // Go to next child:
    if (s.expand_idx != 0) {
      const auto& pexp = exps[s.expand_idx - 1];
      pexp.ctd0.MaybeRemove(res_.tb.ctd);
      pexp.ctd1.MaybeRemove(res_.tb.ctd);
      pexp.pair.MaybeRemove(res_.tb.s);
      if (pexp.idx1.has_value()) unexpanded_.pop_back();
      energy -= pexp.delta;
    }

    if (s.expand_idx != static_cast<int>(exps.size()) && exps[s.expand_idx].delta + energy > delta)
      next_seen = std::min(next_seen, exps[s.expand_idx].delta + energy);

    if (s.expand_idx == static_cast<int>(exps.size()) ||
        exps[s.expand_idx].delta + energy > delta) {
      if (s.should_unexpand) unexpanded_.push_back(to_expand);
      q_.pop_back();
      continue;
    }

    const auto& exp = exps[s.expand_idx++];
    Node ns = {.expand_idx = 0, .to_expand = exp.idx0, .should_unexpand = false};

    // Update global state with this expansion. We can do the others after since
    // they are guaranteed to be empty if this is a terminal.
    energy += exp.delta;
    exp.pair.MaybeApply(res_.tb.s);

    if (!exp.idx0.has_value()) {
      assert(!exp.idx1.has_value());
      assert(!exp.ctd0.IsValid() && !exp.ctd1.IsValid());

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
      if (exp.idx1.has_value()) unexpanded_.push_back(*exp.idx1);
    }

    q_.push_back(ns);
  }
  assert(unexpanded_.empty() && energy == ZERO_E && res_.tb.s == Secondary(res_.tb.s.size()) &&
      res_.tb.ctd == Ctds(res_.tb.ctd.size()));
  return {count, next_seen};
}

template <bool UseLru>
std::vector<Expansion> SuboptIterative<UseLru>::GenerateExpansions(
    const DpIndex& to_expand, Energy delta) const {
  if (std::holds_alternative<base::DpIndex>(to_expand)) {
    auto idx = std::get<base::DpIndex>(to_expand);
    int st = idx.st;
    int en = idx.en;
    int a = idx.a;
    if (en == -1) return ExtExpansions(st, a, delta);
    if (a == DP_P) return PairedOrNoStackExpansions(st, en, /*is_nostack=*/false, delta);
    return UnpairedExpansions(st, en, a, delta);
  }

  if (std::holds_alternative<PenultimateIndex>(to_expand)) {
    auto idx = std::get<PenultimateIndex>(to_expand);
    const int st = idx.st;
    const int en = idx.en;
    const int length = idx.len;
    return PenultimateExpansions(st, en, length, delta);
  }

  auto idx = std::get<NoStackIndex>(to_expand);
  const int st = idx.st;
  const int en = idx.en;
  return PairedOrNoStackExpansions(st, en, /*is_nostack=*/true, delta);
}

template <bool UseLru>
std::vector<Expansion> SuboptIterative<UseLru>::ExtExpansions(int st, int a, Energy delta) const {
  const int N = static_cast<int>(r_.size());
  const auto& dp = dp_.base.dp;
  const auto& ext = dp_.base.ext;
  std::vector<Expansion> exps;
  Energy energy = ZERO_E;

  // Case: No pair starting here
  if (a == EXT) {
    if (st == N) {
      exps.push_back({.delta = ZERO_E});
    } else {
      energy = ext[st + 1][EXT] + m_->pf.Unpaired(st) - ext[st][EXT];
      if (energy <= delta)
        exps.push_back({.delta = energy, .idx0 = base::DpIndex(st + 1, -1, EXT)});
    }
  }
  for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
    // .   .   .   (   .   .   .   )   <   >
    //           stb  st1b   en1b  enb   rem
    const auto stb = r_[st];
    const auto st1b = r_[st + 1];
    const auto enb = r_[en];
    const auto en1b = r_[en - 1];
    const auto base00 = dp[st][en][DP_P] + m_->AuGuPenalty(stb, enb) - ext[st][a];
    const auto base01 = dp[st][en - 1][DP_P] + m_->AuGuPenalty(stb, en1b) - ext[st][a];
    const auto base10 = dp[st + 1][en][DP_P] + m_->AuGuPenalty(st1b, enb) - ext[st][a];
    const auto base11 = dp[st + 1][en - 1][DP_P] + m_->AuGuPenalty(st1b, en1b) - ext[st][a];

    // (   )<.( * ). > Right coax backward
    if (a == EXT_RC && m_->cfg().UseCoaxialStacking()) {
      // Don't set CTDs here since they will have already been set.
      energy = base11 + m_->MismatchCoaxial(en1b, enb, stb, st1b) + m_->pf.Unpaired(st) +
          m_->pf.Unpaired(en) + ext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = base::DpIndex(en + 1, -1, EXT)});
    }

    if (a == EXT_RC) continue;

    // (   )<   >
    energy = base00 + ext[en + 1][EXT];
    if (energy <= delta) {
      // EXT_WC and EXT_GU will have already had their ctds set.
      Expansion exp{.delta = energy,
          .idx0 = base::DpIndex(st, en, DP_P),
          .idx1 = base::DpIndex(en + 1, -1, EXT)};
      if ((a == EXT_WC && IsWcPair(stb, enb)) || (a == EXT_GU && IsGuPair(stb, enb)))
        exps.push_back(exp);

      if (a == EXT) {
        exp.ctd0 = {st, CTD_UNUSED};
        exps.push_back(exp);
      }
    }

    // Only look at EXT from here on.
    if (a != EXT) continue;

    if (m_->cfg().UseDangleMismatch()) {
      // (   )3<   > 3'
      energy = base01 + m_->dangle3[en1b][enb][stb] + m_->pf.Unpaired(en) + ext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({
            .delta = energy,
            .idx0 = base::DpIndex(st, en - 1, DP_P),
            .idx1 = base::DpIndex(en + 1, -1, EXT),
            .ctd0{st, CTD_3_DANGLE},
        });

      // 5(   )<   > 5'
      energy = base10 + m_->dangle5[enb][stb][st1b] + m_->pf.Unpaired(st) + ext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st + 1, en, DP_P),
            .idx1 = base::DpIndex(en + 1, -1, EXT),
            .ctd0{st + 1, CTD_5_DANGLE}});

      // .(   ).<   > Terminal mismatch
      energy = base11 + m_->terminal[en1b][enb][stb][st1b] + m_->pf.Unpaired(st) +
          m_->pf.Unpaired(en) + ext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = base::DpIndex(en + 1, -1, EXT),
            .ctd0{st + 1, CTD_MISMATCH}});
    }

    if (en < N - 1 && m_->cfg().UseCoaxialStacking()) {
      // .(   ).<(   ) > Left coax
      energy = base11 + m_->MismatchCoaxial(en1b, enb, stb, st1b) + m_->pf.Unpaired(st) +
          m_->pf.Unpaired(en);
      if (energy + ext[en + 1][EXT_WC] <= delta)
        exps.push_back({.delta = energy + ext[en + 1][EXT_WC],
            .idx0 = base::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = base::DpIndex(en + 1, -1, EXT_WC),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{en + 1, CTD_LCOAX_WITH_PREV}});

      if (energy + ext[en + 1][EXT_GU] <= delta)
        exps.push_back({.delta = energy + ext[en + 1][EXT_GU],
            .idx0 = base::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = base::DpIndex(en + 1, -1, EXT_GU),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{en + 1, CTD_LCOAX_WITH_PREV}});

      // (   )<.(   ). > Right coax forward
      if (en < N - 2 && base00 + ext[en + 1][EXT_RC] <= delta)
        exps.push_back({.delta = base00 + ext[en + 1][EXT_RC],
            .idx0 = base::DpIndex(st, en, DP_P),
            .idx1 = base::DpIndex(en + 1, -1, EXT_RC),
            .ctd0{st, CTD_RC_WITH_NEXT},
            .ctd1{en + 2, CTD_RC_WITH_PREV}});

      // (   )(<   ) > Flush coax
      energy = base01 + m_->stack[en1b][enb][WcPair(enb)][stb] + ext[en][EXT_WC];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st, en - 1, DP_P),
            .idx1 = base::DpIndex(en, -1, EXT_WC),
            .ctd0{st, CTD_FCOAX_WITH_NEXT},
            .ctd1{en, CTD_FCOAX_WITH_PREV}});

      if (IsGu(enb)) {
        energy = base01 + m_->stack[en1b][enb][GuPair(enb)][stb] + ext[en][EXT_GU];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0 = base::DpIndex(st, en - 1, DP_P),
              .idx1 = base::DpIndex(en, -1, EXT_GU),
              .ctd0{st, CTD_FCOAX_WITH_NEXT},
              .ctd1{en, CTD_FCOAX_WITH_PREV}});
      }
    }
  }

  return exps;
}

template <bool UseLru>
std::vector<Expansion> SuboptIterative<UseLru>::PairedOrNoStackExpansions(
    int st, int en, bool is_nostack, Energy delta) const {
  const auto& dp = dp_.base.dp;
  const auto& nostack = dp_.nostack;
  const auto& penult = dp_.penult;
  std::vector<Expansion> exps;
  Energy energy = ZERO_E;

  const auto stb = r_[st];
  const auto st1b = r_[st + 1];
  const auto st2b = r_[st + 2];
  const auto enb = r_[en];
  const auto en1b = r_[en - 1];
  const auto en2b = r_[en - 2];

  if (!is_nostack) {
    const int max_stack = en - st - HAIRPIN_MIN_SZ + 1;
    const Energy bulge_left = m_->Bulge(r_, st, en, st + 2, en - 1);
    const Energy bulge_right = m_->Bulge(r_, st, en, st + 1, en - 2);

    const auto none = m_->stack[r_[st]][r_[st + 1]][r_[en - 1]][r_[en]] +
        m_->penultimate_stack[en1b][enb][stb][st1b] + m_->pf.Paired(st, en) - dp[st][en][DP_P];
    const auto left = bulge_left + m_->penultimate_stack[en1b][enb][stb][st2b] - dp[st][en][DP_P];
    const auto right = bulge_right + m_->penultimate_stack[en2b][enb][stb][st1b] - dp[st][en][DP_P];

    for (int length = 2; 2 * length <= max_stack; ++length) {
      if (length == 2) {
        energy = none + nostack[st + 1][en - 1] +
            m_->penultimate_stack[r_[st]][r_[st + 1]][r_[en - 1]][r_[en]];
        if (energy <= delta)
          exps.push_back({.delta = energy, .idx0 = NoStackIndex(st + 1, en - 1), .pair{st, en}});
      }

      energy = none + penult[st + 1][en - 1][length - 1];
      if (energy <= delta)
        exps.push_back(
            {.delta = energy, .idx0 = PenultimateIndex(st + 1, en - 1, length - 1), .pair{st, en}});

      if (length == 2) {
        energy = left + nostack[st + 2][en - 1] +
            m_->penultimate_stack[r_[st]][r_[st + 2]][r_[en - 1]][r_[en]];
        if (energy <= delta)
          exps.push_back({.delta = energy, .idx0 = NoStackIndex(st + 2, en - 1), .pair{st, en}});
      }

      energy = left + penult[st + 2][en - 1][length - 1];
      if (energy <= delta)
        exps.push_back(
            {.delta = energy, .idx0 = PenultimateIndex(st + 2, en - 1, length - 1), .pair{st, en}});

      if (length == 2) {
        energy = right + nostack[st + 1][en - 2] +
            m_->penultimate_stack[r_[st]][r_[st + 1]][r_[en - 2]][r_[en]];
        if (energy <= delta)
          exps.push_back({.delta = energy, .idx0 = NoStackIndex(st + 1, en - 2), .pair{st, en}});
      }

      energy = right + penult[st + 1][en - 2][length - 1];
      if (energy <= delta)
        exps.push_back(
            {.delta = energy, .idx0 = PenultimateIndex(st + 1, en - 2, length - 1), .pair{st, en}});
    }
  }

  const auto target = is_nostack ? nostack[st][en] : dp[st][en][DP_P];

  const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
  for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
    for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
      // Try all internal loops. We don't check stacks or 1 nuc bulge loops.
      if (dp[ist][ien][DP_P] < CAP_E && ist - st + en - ien > 3) {
        energy = m_->TwoLoop(r_, st, en, ist, ien) + dp[ist][ien][DP_P] - target;
        if (energy <= delta)
          exps.push_back({.delta = energy, .idx0 = base::DpIndex(ist, ien, DP_P), .pair{st, en}});
      }
    }
  }

  energy = m_->Hairpin(r_, st, en) - target;
  if (energy <= delta) exps.push_back({.delta = energy, .pair{st, en}});

  const auto base_branch_cost = m_->AuGuPenalty(stb, enb) + m_->pf.Paired(st, en) +
      m_->multiloop_a + m_->multiloop_b - target;
  // (<   ><    >)
  energy = base_branch_cost + dp[st + 1][en - 1][DP_U2];
  if (energy <= delta)
    exps.push_back({
        .delta = energy,
        .idx0 = base::DpIndex(st + 1, en - 1, DP_U2),
        .ctd0{en, CTD_UNUSED},
        .pair{st, en},
    });

  if (m_->cfg().UseDangleMismatch()) {
    // (3<   ><   >) 3'
    energy = base_branch_cost + dp[st + 2][en - 1][DP_U2] + m_->dangle3[stb][st1b][enb] +
        m_->pf.Unpaired(st + 1);
    if (energy <= delta)
      exps.push_back({.delta = energy,
          .idx0 = base::DpIndex(st + 2, en - 1, DP_U2),
          .ctd0{en, CTD_3_DANGLE},
          .pair{st, en}});

    // (<   ><   >5) 5'
    energy = base_branch_cost + dp[st + 1][en - 2][DP_U2] + m_->dangle5[stb][en1b][enb] +
        m_->pf.Unpaired(en - 1);
    if (energy <= delta)
      exps.push_back({.delta = energy,
          .idx0 = base::DpIndex(st + 1, en - 2, DP_U2),
          .ctd0{en, CTD_5_DANGLE},
          .pair{st, en}});

    // (.<   ><   >.) Terminal mismatch
    energy = base_branch_cost + dp[st + 2][en - 2][DP_U2] + m_->terminal[stb][st1b][en1b][enb] +
        m_->pf.Unpaired(st + 1) + m_->pf.Unpaired(en - 1);
    if (energy <= delta)
      exps.push_back({.delta = energy,
          .idx0 = base::DpIndex(st + 2, en - 2, DP_U2),
          .ctd0{en, CTD_MISMATCH},
          .pair{st, en}});
  }

  if (m_->cfg().UseCoaxialStacking()) {
    const auto outer_coax = m_->MismatchCoaxial(stb, st1b, en1b, enb) + m_->pf.Unpaired(st + 1) +
        m_->pf.Unpaired(en - 1);
    for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
      const Base pl1b = r_[piv - 1];
      const Base plb = r_[piv];
      const Base prb = r_[piv + 1];
      const Base pr1b = r_[piv + 2];

      // (.(   )   .) Left outer coax - P
      energy = base_branch_cost + dp[st + 2][piv][DP_P] + m_->multiloop_b +
          m_->AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax;
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st + 2, piv, DP_P),
            .idx1 = base::DpIndex(piv + 1, en - 2, DP_U),
            .ctd0{en, CTD_LCOAX_WITH_NEXT},
            .ctd1{st + 2, CTD_LCOAX_WITH_PREV},
            .pair{st, en}});

      // (.   (   ).) Right outer coax
      energy = base_branch_cost + dp[st + 2][piv][DP_U] + m_->multiloop_b +
          m_->AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax;
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st + 2, piv, DP_U),
            .idx1 = base::DpIndex(piv + 1, en - 2, DP_P),
            .ctd0{en, CTD_RC_WITH_PREV},
            .ctd1{piv + 1, CTD_RC_WITH_NEXT},
            .pair{st, en}});

      // (.(   ).   ) Left inner coax
      energy = base_branch_cost + dp[st + 2][piv - 1][DP_P] + m_->multiloop_b +
          m_->AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
          m_->MismatchCoaxial(pl1b, plb, st1b, st2b) + m_->pf.Unpaired(st + 1) +
          m_->pf.Unpaired(piv);
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st + 2, piv - 1, DP_P),
            .idx1 = base::DpIndex(piv + 1, en - 1, DP_U),
            .ctd0{en, CTD_RC_WITH_NEXT},
            .ctd1{st + 2, CTD_RC_WITH_PREV},
            .pair{st, en}});

      // (   .(   ).) Right inner coax
      energy = base_branch_cost + dp[st + 1][piv][DP_U] + m_->multiloop_b +
          m_->AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
          m_->MismatchCoaxial(en2b, en1b, prb, pr1b) + m_->pf.Unpaired(piv + 1) +
          m_->pf.Unpaired(en - 1);
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st + 1, piv, DP_U),
            .idx1 = base::DpIndex(piv + 2, en - 2, DP_P),
            .ctd0{en, CTD_LCOAX_WITH_PREV},
            .ctd1{piv + 2, CTD_LCOAX_WITH_NEXT},
            .pair{st, en}});

      // ((   )   ) Left flush coax
      energy = base_branch_cost + dp[st + 1][piv][DP_P] + m_->multiloop_b +
          m_->AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] + m_->stack[stb][st1b][plb][enb];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st + 1, piv, DP_P),
            .idx1 = base::DpIndex(piv + 1, en - 1, DP_U),
            .ctd0{en, CTD_FCOAX_WITH_NEXT},
            .ctd1{st + 1, CTD_FCOAX_WITH_PREV},
            .pair{st, en}});

      // (   (   )) Right flush coax
      energy = base_branch_cost + dp[st + 1][piv][DP_U] + m_->multiloop_b +
          m_->AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] + m_->stack[stb][prb][en1b][enb];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st + 1, piv, DP_U),
            .idx1 = base::DpIndex(piv + 1, en - 1, DP_P),
            .ctd0{en, CTD_FCOAX_WITH_PREV},
            .ctd1{piv + 1, CTD_FCOAX_WITH_NEXT},
            .pair{st, en}});
    }
  }

  return exps;
}

template <bool UseLru>
std::vector<Expansion> SuboptIterative<UseLru>::UnpairedExpansions(
    int st, int en, int a, Energy delta) const {
  const auto& dp = dp_.base.dp;
  std::vector<Expansion> exps;
  Energy energy = ZERO_E;

  const auto stb = r_[st];
  const auto st1b = r_[st + 1];

  // Left unpaired. Either DP_U or DP_U2.
  if (st + 1 < en && (a == DP_U || a == DP_U2)) {
    energy = dp[st + 1][en][a] + m_->pf.Unpaired(st) - dp[st][en][a];
    if (energy <= delta) exps.push_back({.delta = energy, .idx0 = base::DpIndex(st + 1, en, a)});
  }

  // Pair here.
  for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
    //   (   .   )<   (
    // stb pl1b pb   pr1b
    const auto pb = r_[piv];
    const auto pl1b = r_[piv - 1];
    // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the
    // right.
    const auto base00 =
        dp[st][piv][DP_P] + m_->AuGuPenalty(stb, pb) + m_->multiloop_b - dp[st][en][a];
    const auto base01 =
        dp[st][piv - 1][DP_P] + m_->AuGuPenalty(stb, pl1b) + m_->multiloop_b - dp[st][en][a];
    const auto base10 =
        dp[st + 1][piv][DP_P] + m_->AuGuPenalty(st1b, pb) + m_->multiloop_b - dp[st][en][a];
    const auto base11 =
        dp[st + 1][piv - 1][DP_P] + m_->AuGuPenalty(st1b, pl1b) + m_->multiloop_b - dp[st][en][a];

    const auto right_paired = dp[piv + 1][en][DP_U];
    // This is only usable if a != DP_U2 since this leaves everything unpaired.
    const auto right_unpaired = m_->pf.UnpairedCum(piv + 1, en);

    if (m_->cfg().UseCoaxialStacking()) {
      // Check a == U_RC:
      // (   )<.( ** ). > Right coax backward
      if (a == DP_U_RC) {
        energy = base11 + m_->MismatchCoaxial(pl1b, pb, stb, st1b) + m_->pf.Unpaired(st) +
            m_->pf.Unpaired(piv);
        if (energy + right_unpaired <= delta)
          exps.push_back(
              {.delta = energy + right_unpaired, .idx0 = base::DpIndex(st + 1, piv - 1, DP_P)});

        if (energy + right_paired <= delta)
          exps.push_back({.delta = energy + right_paired,
              .idx0 = base::DpIndex(st + 1, piv - 1, DP_P),
              .idx1 = base::DpIndex(piv + 1, en, DP_U)});
      }
    }

    // DP_U_RC is only the above case.
    if (a == DP_U_RC) continue;

    // (   )<   > - U, U2, U_WC?, U_GU?
    if (a == DP_U) {
      energy = base00 + right_unpaired;
      if (energy <= delta)
        exps.push_back(
            {.delta = energy, .idx0 = base::DpIndex(st, piv, DP_P), .ctd0 = {st, CTD_UNUSED}});

      energy = base00 + right_paired;
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st, piv, DP_P),
            .idx1 = base::DpIndex(piv + 1, en, DP_U),
            .ctd0 = {st, CTD_UNUSED}});
    }

    energy = base00 + right_paired;
    if (a == DP_U2 && energy <= delta)
      exps.push_back({.delta = energy,
          .idx0 = base::DpIndex(st, piv, DP_P),
          .idx1 = base::DpIndex(piv + 1, en, DP_U),
          .ctd0 = {st, CTD_UNUSED}});

    if ((a == DP_U_WC && IsWcPair(stb, pb)) || (a == DP_U_GU && IsGuPair(stb, pb))) {
      energy = base00 + right_unpaired;
      if (energy <= delta) exps.push_back({.delta = energy, .idx0 = base::DpIndex(st, piv, DP_P)});

      energy = base00 + right_paired;
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = base::DpIndex(st, piv, DP_P),
            .idx1 = base::DpIndex(piv + 1, en, DP_U)});
    }

    // The rest of the cases are for U and U2.
    if (a != DP_U && a != DP_U2) continue;

    if (m_->cfg().UseDangleMismatch()) {
      // (   )3<   > 3' - U, U2
      energy = base01 + m_->dangle3[pl1b][pb][stb] + m_->pf.Unpaired(piv);
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired,
            .idx0 = base::DpIndex(st, piv - 1, DP_P),
            .ctd0{st, CTD_3_DANGLE}});

      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0 = base::DpIndex(st, piv - 1, DP_P),
            .idx1 = base::DpIndex(piv + 1, en, DP_U),
            .ctd0{st, CTD_3_DANGLE}});

      // 5(   )<   > 5' - U, U2
      energy = base10 + m_->dangle5[pb][stb][st1b] + m_->pf.Unpaired(st);
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired,
            .idx0 = base::DpIndex(st + 1, piv, DP_P),
            .ctd0{st + 1, CTD_5_DANGLE}});

      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0 = base::DpIndex(st + 1, piv, DP_P),
            .idx1 = base::DpIndex(piv + 1, en, DP_U),
            .ctd0{st + 1, CTD_5_DANGLE}});

      // .(   ).<   > Terminal mismatch - U, U2
      energy =
          base11 + m_->terminal[pl1b][pb][stb][st1b] + m_->pf.Unpaired(st) + m_->pf.Unpaired(piv);
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired,
            .idx0 = base::DpIndex(st + 1, piv - 1, DP_P),
            .ctd0{st + 1, CTD_MISMATCH}});

      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0 = base::DpIndex(st + 1, piv - 1, DP_P),
            .idx1 = base::DpIndex(piv + 1, en, DP_U),
            .ctd0{st + 1, CTD_MISMATCH}});
    }

    if (m_->cfg().UseCoaxialStacking()) {
      // .(   ).<(   ) > Left coax - U, U2
      energy = base11 + m_->MismatchCoaxial(pl1b, pb, stb, st1b) + m_->pf.Unpaired(st) +
          m_->pf.Unpaired(piv);
      if (energy + dp[piv + 1][en][DP_U_WC] <= delta)
        exps.push_back({
            .delta = energy + dp[piv + 1][en][DP_U_WC],
            .idx0 = base::DpIndex(st + 1, piv - 1, DP_P),
            .idx1 = base::DpIndex(piv + 1, en, DP_U_WC),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{piv + 1, CTD_LCOAX_WITH_PREV},
        });
      if (energy + dp[piv + 1][en][DP_U_GU] <= delta)
        exps.push_back({
            .delta = energy + dp[piv + 1][en][DP_U_GU],
            .idx0 = base::DpIndex(st + 1, piv - 1, DP_P),
            .idx1 = base::DpIndex(piv + 1, en, DP_U_GU),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{piv + 1, CTD_LCOAX_WITH_PREV},
        });

      // (   )<.(   ). > Right coax forward - U, U2
      energy = base00 + dp[piv + 1][en][DP_U_RC];
      if (energy <= delta)
        exps.push_back({
            .delta = energy,
            .idx0 = base::DpIndex(st, piv, DP_P),
            .idx1 = base::DpIndex(piv + 1, en, DP_U_RC),
            .ctd0{st, CTD_RC_WITH_NEXT},
            .ctd1{piv + 2, CTD_RC_WITH_PREV},
        });

      // (   )(<   ) > Flush coax - U, U2
      energy = base01 + m_->stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC];
      if (energy <= delta)
        exps.push_back({
            .delta = energy,
            .idx0 = base::DpIndex(st, piv - 1, DP_P),
            .idx1 = base::DpIndex(piv, en, DP_U_WC),
            .ctd0{st, CTD_FCOAX_WITH_NEXT},
            .ctd1{piv, CTD_FCOAX_WITH_PREV},
        });

      if ((IsGu(pb))) {
        energy = base01 + m_->stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU];
        if (energy <= delta)
          exps.push_back({
              .delta = energy,
              .idx0 = base::DpIndex(st, piv - 1, DP_P),
              .idx1 = base::DpIndex(piv, en, DP_U_GU),
              .ctd0{st, CTD_FCOAX_WITH_NEXT},
              .ctd1{piv, CTD_FCOAX_WITH_PREV},
          });
      }
    }
  }

  return exps;
}

template <bool UseLru>
std::vector<Expansion> SuboptIterative<UseLru>::PenultimateExpansions(
    int st, int en, int length, Energy delta) const {
  const auto& nostack = dp_.nostack;
  const auto& penult = dp_.penult;
  std::vector<Expansion> exps;
  Energy energy = ZERO_E;

  const auto bulge_left = m_->Bulge(r_, st, en, st + 2, en - 1);
  const auto bulge_right = m_->Bulge(r_, st, en, st + 1, en - 2);

  auto none = m_->stack[r_[st]][r_[st + 1]][r_[en - 1]][r_[en]] + m_->pf.Paired(st, en) -
      penult[st][en][length];
  if (length == 2) {
    energy = none + nostack[st + 1][en - 1] +
        m_->penultimate_stack[r_[st]][r_[st + 1]][r_[en - 1]][r_[en]];
    if (energy <= delta)
      exps.push_back({.delta = energy, .idx0 = NoStackIndex(st + 1, en - 1), .pair{st, en}});
  }
  energy = none + penult[st + 1][en - 1][length - 1];
  if (energy <= delta)
    exps.push_back(
        {.delta = energy, .idx0 = PenultimateIndex(st + 1, en - 1, length - 1), .pair{st, en}});

  auto left = bulge_left - penult[st][en][length];
  if (length == 2) {
    energy = left + nostack[st + 2][en - 1] +
        m_->penultimate_stack[r_[st]][r_[st + 2]][r_[en - 1]][r_[en]];
    if (energy <= delta)
      exps.push_back({.delta = energy, .idx0 = NoStackIndex(st + 2, en - 1), .pair{st, en}});
  }

  energy = left + penult[st + 2][en - 1][length - 1];
  if (energy <= delta)
    exps.push_back(
        {.delta = energy, .idx0 = PenultimateIndex(st + 2, en - 1, length - 1), .pair{st, en}});

  auto right = bulge_right - penult[st][en][length];
  if (length == 2) {
    energy = right + nostack[st + 1][en - 2] +
        m_->penultimate_stack[r_[st]][r_[st + 1]][r_[en - 2]][r_[en]];
    if (energy <= delta)
      exps.push_back({.delta = energy, .idx0 = NoStackIndex(st + 1, en - 2), .pair{st, en}});
  }

  energy = right + penult[st + 1][en - 2][length - 1];
  if (energy <= delta)
    exps.push_back(
        {.delta = energy, .idx0 = PenultimateIndex(st + 1, en - 2, length - 1), .pair{st, en}});

  return exps;
}

template class SuboptIterative<false>;
template class SuboptIterative<true>;

}  // namespace mrna::md::stack
