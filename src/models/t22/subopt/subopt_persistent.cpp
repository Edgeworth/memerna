// Copyright 2023 E.
#include "models/t22/subopt/subopt_persistent.h"

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
using t04::DP_U;
using t04::DP_U2;
using t04::DP_U_GU;
using t04::DP_U_RC;
using t04::DP_U_WC;
using t04::EXT;
using t04::EXT_GU;
using t04::EXT_RC;
using t04::EXT_WC;

namespace {

// Check the time every `CHECK_TIME_FREQ` iterations.
constexpr int CHECK_TIME_FREQ = 10000;

}  // namespace

SuboptPersistent::SuboptPersistent(Primary r, Model::Ptr em, DpState dp, SuboptCfg cfg)
    : r_(std::move(r)), em_(std::move(em)), dp_(std::move(dp)), cfg_(cfg) {}

int SuboptPersistent::Run(const SuboptCallback& fn) {
  res_ = SuboptResult(ZERO_E, trace::TraceResult(Secondary(r_.size()), Ctds(r_.size())));
  q_.reserve(r_.size());
  cache_.resize(MaxLinearIndex(r_.size()));

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false, true},
      .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::NONE},
  };
  support.VerifySupported(funcname(), em_->cfg);

  spdlog::debug("t22 {} with cfg {}", funcname(), em_->cfg);

  q_.clear();
  pq_ = {};  // priority queue has no clear method
  q_.reserve(r_.size() * r_.size());
  cache_.resize(MaxLinearIndex(r_.size()));

  const auto start_idx = t04::DpIndex(0, -1, EXT);
  const Energy mfe = dp_.t04.Index(start_idx);
  q_.push_back({.expand_idx = 0, .to_expand = start_idx});
  pq_.emplace(0, 0);

  int num_strucs = 0;
  auto start_time = std::chrono::steady_clock::now();

  while (!pq_.empty()) {
    if (num_strucs >= cfg_.strucs) break;

    if (cfg_.time_secs >= 0.0 && num_strucs % CHECK_TIME_FREQ == 0) {
      auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - start_time);
      if (elapsed.count() >= cfg_.time_secs) break;
    }
    auto [delta, idx] = RunInternal();
    if (idx == -1 || delta > cfg_.delta) break;

    num_strucs++;
    GenerateResult(idx);
    res_.energy = delta + mfe;
    fn(res_);
  }

  return num_strucs;
}

std::pair<Energy, int> SuboptPersistent::RunInternal() {
  while (!pq_.empty()) {
    auto [neg_delta, idx] = pq_.top();
    pq_.pop();
    auto& s = q_[idx];

    if (!s.to_expand.has_value()) {
      return {-neg_delta, idx};
    }

    const auto& exps = GetExpansion(*s.to_expand);
    assert(!exps.empty());  // Must produce at least one expansion: {-1, -1, -1}.

    const auto& exp = exps[s.expand_idx];
    const bool has_unexpanded = exp.idx1.has_value();
    Node ns = {.parent_idx = idx,
        .parent_expand_idx = s.expand_idx,
        .unexpanded_idx = has_unexpanded ? idx : s.unexpanded_idx,
        .unexpanded_expand_idx = has_unexpanded ? s.expand_idx : s.unexpanded_expand_idx,
        .expand_idx = 0,
        .to_expand = exp.idx0};
    s.expand_idx++;

    if (s.expand_idx != static_cast<int>(exps.size())) {
      pq_.emplace(neg_delta + exp.delta - exps[s.expand_idx].delta, idx);
    }

    if (!exp.idx0.has_value() && s.unexpanded_idx != -1) {
      assert(!has_unexpanded);
      assert(!exp.ctd0.IsValid() && !exp.ctd1.IsValid());
      assert(q_[s.unexpanded_idx].to_expand.has_value());

      const auto& unexpanded_exp =
          GetExpansion(*q_[s.unexpanded_idx].to_expand)[s.unexpanded_expand_idx];
      assert(unexpanded_exp.idx1.has_value());  // There should be a valid unexpanded DpIndex here.

      ns.unexpanded_expand_idx = q_[ns.unexpanded_idx].unexpanded_expand_idx;
      ns.unexpanded_idx = q_[ns.unexpanded_idx].unexpanded_idx;
      ns.to_expand = unexpanded_exp.idx1;
    }

    pq_.emplace(neg_delta, static_cast<int>(q_.size()));
    q_.push_back(ns);
  }
  return {ZERO_E, -1};
}

void SuboptPersistent::GenerateResult(int idx) {
  res_.tb.s.reset(r_.size());
  res_.tb.ctd.reset(r_.size());

  int expand_idx = q_[idx].parent_expand_idx;
  idx = q_[idx].parent_idx;
  while (idx != -1) {
    const auto to_expand = q_[idx].to_expand;
    assert(to_expand.has_value());
    const auto& exp = GetExpansion(*to_expand)[expand_idx];

    exp.pair.MaybeApply(res_.tb.s);
    exp.ctd0.MaybeApply(res_.tb.ctd);
    exp.ctd1.MaybeApply(res_.tb.ctd);

    expand_idx = q_[idx].parent_expand_idx;
    idx = q_[idx].parent_idx;
  }
}

std::vector<Expansion> SuboptPersistent::GenerateExpansions(
    const DpIndex& to_expand, Energy delta) const {
  if (std::holds_alternative<t04::DpIndex>(to_expand)) {
    auto idx = std::get<t04::DpIndex>(to_expand);
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

std::vector<Expansion> SuboptPersistent::ExtExpansions(int st, int a, Energy delta) const {
  const int N = static_cast<int>(r_.size());
  const auto& dp = dp_.t04.dp;
  const auto& ext = dp_.t04.ext;
  std::vector<Expansion> exps;
  Energy energy = ZERO_E;

  // Case: No pair starting here
  if (a == EXT) {
    if (st == N) {
      exps.push_back({.delta = ZERO_E});
    } else {
      energy = ext[st + 1][EXT] + em_->PfUnpaired(st) - ext[st][EXT];
      if (energy <= delta) exps.push_back({.delta = energy, .idx0 = t04::DpIndex(st + 1, -1, EXT)});
    }
  }
  for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
    // .   .   .   (   .   .   .   )   <   >
    //           stb  st1b   en1b  enb   rem
    const auto stb = r_[st];
    const auto st1b = r_[st + 1];
    const auto enb = r_[en];
    const auto en1b = r_[en - 1];
    const auto base00 = dp[st][en][DP_P] + em_->AuGuPenalty(stb, enb) - ext[st][a];
    const auto base01 = dp[st][en - 1][DP_P] + em_->AuGuPenalty(stb, en1b) - ext[st][a];
    const auto base10 = dp[st + 1][en][DP_P] + em_->AuGuPenalty(st1b, enb) - ext[st][a];
    const auto base11 = dp[st + 1][en - 1][DP_P] + em_->AuGuPenalty(st1b, en1b) - ext[st][a];

    // (   )<.( * ). > Right coax backward
    if (a == EXT_RC && em_->cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
      // Don't set CTDs here since they will have already been set.
      energy = base11 + em_->MismatchCoaxial(en1b, enb, stb, st1b) + em_->PfUnpaired(st) +
          em_->PfUnpaired(en) + ext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT)});
    }

    if (a == EXT_RC) continue;

    // (   )<   >
    energy = base00 + ext[en + 1][EXT];
    if (energy <= delta) {
      // EXT_WC and EXT_GU will have already had their ctds set.
      Expansion exp{.delta = energy,
          .idx0 = t04::DpIndex(st, en, DP_P),
          .idx1 = t04::DpIndex(en + 1, -1, EXT)};
      if ((a == EXT_WC && IsWcPair(stb, enb)) || (a == EXT_GU && IsGuPair(stb, enb)))
        exps.push_back(exp);

      if (a == EXT) {
        exp.ctd0 = {st, CTD_UNUSED};
        exps.push_back(exp);
      }
    }

    // Only look at EXT from here on.
    if (a != EXT) continue;

    if (em_->cfg.ctd == erg::EnergyCfg::Ctd::ALL || em_->cfg.ctd == erg::EnergyCfg::Ctd::NO_COAX) {
      // (   )3<   > 3'
      energy = base01 + em_->dangle3[en1b][enb][stb] + em_->PfUnpaired(en) + ext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({
            .delta = energy,
            .idx0 = t04::DpIndex(st, en - 1, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT),
            .ctd0{st, CTD_3_DANGLE},
        });

      // 5(   )<   > 5'
      energy = base10 + em_->dangle5[enb][stb][st1b] + em_->PfUnpaired(st) + ext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st + 1, en, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT),
            .ctd0{st + 1, CTD_5_DANGLE}});

      // .(   ).<   > Terminal mismatch
      energy = base11 + em_->terminal[en1b][enb][stb][st1b] + em_->PfUnpaired(st) +
          em_->PfUnpaired(en) + ext[en + 1][EXT];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT),
            .ctd0{st + 1, CTD_MISMATCH}});
    }

    if (en < N - 1 && em_->cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
      // .(   ).<(   ) > Left coax
      energy = base11 + em_->MismatchCoaxial(en1b, enb, stb, st1b) + em_->PfUnpaired(st) +
          em_->PfUnpaired(en);
      if (energy + ext[en + 1][EXT_WC] <= delta)
        exps.push_back({.delta = energy + ext[en + 1][EXT_WC],
            .idx0 = t04::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT_WC),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{en + 1, CTD_LCOAX_WITH_PREV}});

      if (energy + ext[en + 1][EXT_GU] <= delta)
        exps.push_back({.delta = energy + ext[en + 1][EXT_GU],
            .idx0 = t04::DpIndex(st + 1, en - 1, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT_GU),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{en + 1, CTD_LCOAX_WITH_PREV}});

      // (   )<.(   ). > Right coax forward
      if (en < N - 2 && base00 + ext[en + 1][EXT_RC] <= delta)
        exps.push_back({.delta = base00 + ext[en + 1][EXT_RC],
            .idx0 = t04::DpIndex(st, en, DP_P),
            .idx1 = t04::DpIndex(en + 1, -1, EXT_RC),
            .ctd0{st, CTD_RC_WITH_NEXT},
            .ctd1{en + 2, CTD_RC_WITH_PREV}});

      // (   )(<   ) > Flush coax
      energy = base01 + em_->stack[en1b][enb][WcPair(enb)][stb] + ext[en][EXT_WC];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st, en - 1, DP_P),
            .idx1 = t04::DpIndex(en, -1, EXT_WC),
            .ctd0{st, CTD_FCOAX_WITH_NEXT},
            .ctd1{en, CTD_FCOAX_WITH_PREV}});

      if (IsGu(enb)) {
        energy = base01 + em_->stack[en1b][enb][GuPair(enb)][stb] + ext[en][EXT_GU];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0 = t04::DpIndex(st, en - 1, DP_P),
              .idx1 = t04::DpIndex(en, -1, EXT_GU),
              .ctd0{st, CTD_FCOAX_WITH_NEXT},
              .ctd1{en, CTD_FCOAX_WITH_PREV}});
      }
    }
  }

  return exps;
}

std::vector<Expansion> SuboptPersistent::PairedOrNoStackExpansions(
    int st, int en, bool is_nostack, Energy delta) const {
  const auto& dp = dp_.t04.dp;
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
    const Energy bulge_left = em_->Bulge(r_, st, en, st + 2, en - 1);
    const Energy bulge_right = em_->Bulge(r_, st, en, st + 1, en - 2);

    const auto none = em_->stack[r_[st]][r_[st + 1]][r_[en - 1]][r_[en]] +
        em_->penultimate_stack[en1b][enb][stb][st1b] + em_->PfPaired(st, en) - dp[st][en][DP_P];
    const auto left = bulge_left + em_->penultimate_stack[en1b][enb][stb][st2b] - dp[st][en][DP_P];
    const auto right =
        bulge_right + em_->penultimate_stack[en2b][enb][stb][st1b] - dp[st][en][DP_P];

    for (int length = 2; 2 * length <= max_stack; ++length) {
      if (length == 2) {
        energy = none + nostack[st + 1][en - 1] +
            em_->penultimate_stack[r_[st]][r_[st + 1]][r_[en - 1]][r_[en]];
        if (energy <= delta)
          exps.push_back({.delta = energy, .idx0 = NoStackIndex(st + 1, en - 1), .pair{st, en}});
      }

      energy = none + penult[st + 1][en - 1][length - 1];
      if (energy <= delta)
        exps.push_back(
            {.delta = energy, .idx0 = PenultimateIndex(st + 1, en - 1, length - 1), .pair{st, en}});

      if (length == 2) {
        energy = left + nostack[st + 2][en - 1] +
            em_->penultimate_stack[r_[st]][r_[st + 2]][r_[en - 1]][r_[en]];
        if (energy <= delta)
          exps.push_back({.delta = energy, .idx0 = NoStackIndex(st + 2, en - 1), .pair{st, en}});
      }

      energy = left + penult[st + 2][en - 1][length - 1];
      if (energy <= delta)
        exps.push_back(
            {.delta = energy, .idx0 = PenultimateIndex(st + 2, en - 1, length - 1), .pair{st, en}});

      if (length == 2) {
        energy = right + nostack[st + 1][en - 2] +
            em_->penultimate_stack[r_[st]][r_[st + 1]][r_[en - 2]][r_[en]];
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
        energy = em_->TwoLoop(r_, st, en, ist, ien) + dp[ist][ien][DP_P] - target;
        if (energy <= delta)
          exps.push_back({.delta = energy, .idx0 = t04::DpIndex(ist, ien, DP_P), .pair{st, en}});
      }
    }
  }

  energy = em_->Hairpin(r_, st, en) - target;
  if (energy <= delta) exps.push_back({.delta = energy, .pair{st, en}});

  const auto base_branch_cost = em_->AuGuPenalty(stb, enb) + em_->PfPaired(st, en) +
      em_->multiloop_hack_a + em_->multiloop_hack_b - target;
  // (<   ><    >)
  energy = base_branch_cost + dp[st + 1][en - 1][DP_U2];
  if (energy <= delta)
    exps.push_back({
        .delta = energy,
        .idx0 = t04::DpIndex(st + 1, en - 1, DP_U2),
        .ctd0{en, CTD_UNUSED},
        .pair{st, en},
    });

  if (em_->cfg.ctd == erg::EnergyCfg::Ctd::ALL || em_->cfg.ctd == erg::EnergyCfg::Ctd::NO_COAX) {
    // (3<   ><   >) 3'
    energy = base_branch_cost + dp[st + 2][en - 1][DP_U2] + em_->dangle3[stb][st1b][enb] +
        em_->PfUnpaired(st + 1);
    if (energy <= delta)
      exps.push_back({.delta = energy,
          .idx0 = t04::DpIndex(st + 2, en - 1, DP_U2),
          .ctd0{en, CTD_3_DANGLE},
          .pair{st, en}});

    // (<   ><   >5) 5'
    energy = base_branch_cost + dp[st + 1][en - 2][DP_U2] + em_->dangle5[stb][en1b][enb] +
        em_->PfUnpaired(en - 1);
    if (energy <= delta)
      exps.push_back({.delta = energy,
          .idx0 = t04::DpIndex(st + 1, en - 2, DP_U2),
          .ctd0{en, CTD_5_DANGLE},
          .pair{st, en}});

    // (.<   ><   >.) Terminal mismatch
    energy = base_branch_cost + dp[st + 2][en - 2][DP_U2] + em_->terminal[stb][st1b][en1b][enb] +
        em_->PfUnpaired(st + 1) + em_->PfUnpaired(en - 1);
    if (energy <= delta)
      exps.push_back({.delta = energy,
          .idx0 = t04::DpIndex(st + 2, en - 2, DP_U2),
          .ctd0{en, CTD_MISMATCH},
          .pair{st, en}});
  }

  if (em_->cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
    for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
      const Base pl1b = r_[piv - 1];
      const Base plb = r_[piv];
      const Base prb = r_[piv + 1];
      const Base pr1b = r_[piv + 2];

      // (.(   )   .) Left outer coax - P
      const auto outer_coax = em_->MismatchCoaxial(stb, st1b, en1b, enb) + em_->PfUnpaired(st + 1) +
          em_->PfUnpaired(en - 1);

      energy = base_branch_cost + dp[st + 2][piv][DP_P] + em_->multiloop_hack_b +
          em_->AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax;
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st + 2, piv, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en - 2, DP_U),
            .ctd0{en, CTD_LCOAX_WITH_NEXT},
            .ctd1{st + 2, CTD_LCOAX_WITH_PREV},
            .pair{st, en}});

      // (.   (   ).) Right outer coax
      energy = base_branch_cost + dp[st + 2][piv][DP_U] + em_->multiloop_hack_b +
          em_->AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax;
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st + 2, piv, DP_U),
            .idx1 = t04::DpIndex(piv + 1, en - 2, DP_P),
            .ctd0{en, CTD_RC_WITH_PREV},
            .ctd1{piv + 1, CTD_RC_WITH_NEXT},
            .pair{st, en}});

      // (.(   ).   ) Left inner coax
      energy = base_branch_cost + dp[st + 2][piv - 1][DP_P] + em_->multiloop_hack_b +
          em_->AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
          em_->MismatchCoaxial(pl1b, plb, st1b, st2b) + em_->PfUnpaired(st + 1) +
          em_->PfUnpaired(piv);
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st + 2, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en - 1, DP_U),
            .ctd0{en, CTD_RC_WITH_NEXT},
            .ctd1{st + 2, CTD_RC_WITH_PREV},
            .pair{st, en}});

      // (   .(   ).) Right inner coax
      energy = base_branch_cost + dp[st + 1][piv][DP_U] + em_->multiloop_hack_b +
          em_->AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
          em_->MismatchCoaxial(en2b, en1b, prb, pr1b) + em_->PfUnpaired(piv + 1) +
          em_->PfUnpaired(en - 1);
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st + 1, piv, DP_U),
            .idx1 = t04::DpIndex(piv + 2, en - 2, DP_P),
            .ctd0{en, CTD_LCOAX_WITH_PREV},
            .ctd1{piv + 2, CTD_LCOAX_WITH_NEXT},
            .pair{st, en}});

      // ((   )   ) Left flush coax
      energy = base_branch_cost + dp[st + 1][piv][DP_P] + em_->multiloop_hack_b +
          em_->AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] + em_->stack[stb][st1b][plb][enb];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st + 1, piv, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en - 1, DP_U),
            .ctd0{en, CTD_FCOAX_WITH_NEXT},
            .ctd1{st + 1, CTD_FCOAX_WITH_PREV},
            .pair{st, en}});

      // (   (   )) Right flush coax
      energy = base_branch_cost + dp[st + 1][piv][DP_U] + em_->multiloop_hack_b +
          em_->AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] + em_->stack[stb][prb][en1b][enb];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st + 1, piv, DP_U),
            .idx1 = t04::DpIndex(piv + 1, en - 1, DP_P),
            .ctd0{en, CTD_FCOAX_WITH_PREV},
            .ctd1{piv + 1, CTD_FCOAX_WITH_NEXT},
            .pair{st, en}});
    }
  }

  return exps;
}

std::vector<Expansion> SuboptPersistent::UnpairedExpansions(
    int st, int en, int a, Energy delta) const {
  const auto& dp = dp_.t04.dp;
  std::vector<Expansion> exps;
  Energy energy = ZERO_E;

  const auto stb = r_[st];
  const auto st1b = r_[st + 1];

  // Left unpaired. Either DP_U or DP_U2.
  if (st + 1 < en && (a == DP_U || a == DP_U2)) {
    energy = dp[st + 1][en][a] + em_->PfUnpaired(st) - dp[st][en][a];
    if (energy <= delta) exps.push_back({.delta = energy, .idx0 = t04::DpIndex(st + 1, en, a)});
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
        dp[st][piv][DP_P] + em_->AuGuPenalty(stb, pb) + em_->multiloop_hack_b - dp[st][en][a];
    const auto base01 =
        dp[st][piv - 1][DP_P] + em_->AuGuPenalty(stb, pl1b) + em_->multiloop_hack_b - dp[st][en][a];
    const auto base10 =
        dp[st + 1][piv][DP_P] + em_->AuGuPenalty(st1b, pb) + em_->multiloop_hack_b - dp[st][en][a];
    const auto base11 = dp[st + 1][piv - 1][DP_P] + em_->AuGuPenalty(st1b, pl1b) +
        em_->multiloop_hack_b - dp[st][en][a];

    const auto right_paired = dp[piv + 1][en][DP_U];
    // This is only usable if a != DP_U2 since this leaves everything unpaired.
    const auto right_unpaired = em_->PfUnpairedCum(piv + 1, en);

    if (em_->cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
      // Check a == U_RC:
      // (   )<.( ** ). > Right coax backward
      if (a == DP_U_RC) {
        energy = base11 + em_->MismatchCoaxial(pl1b, pb, stb, st1b) + em_->PfUnpaired(st) +
            em_->PfUnpaired(piv);
        if (energy + right_unpaired <= delta)
          exps.push_back(
              {.delta = energy + right_unpaired, .idx0 = t04::DpIndex(st + 1, piv - 1, DP_P)});

        if (energy + right_paired <= delta)
          exps.push_back({.delta = energy + right_paired,
              .idx0 = t04::DpIndex(st + 1, piv - 1, DP_P),
              .idx1 = t04::DpIndex(piv + 1, en, DP_U)});
      }
    }

    // DP_U_RC is only the above case.
    if (a == DP_U_RC) continue;

    // (   )<   > - U, U2, U_WC?, U_GU?
    if (a == DP_U) {
      energy = base00 + right_unpaired;
      if (energy <= delta)
        exps.push_back(
            {.delta = energy, .idx0 = t04::DpIndex(st, piv, DP_P), .ctd0 = {st, CTD_UNUSED}});

      energy = base00 + right_paired;
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st, piv, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U),
            .ctd0 = {st, CTD_UNUSED}});
    }

    energy = base00 + right_paired;
    if (a == DP_U2 && energy <= delta)
      exps.push_back({.delta = energy,
          .idx0 = t04::DpIndex(st, piv, DP_P),
          .idx1 = t04::DpIndex(piv + 1, en, DP_U),
          .ctd0 = {st, CTD_UNUSED}});

    if ((a == DP_U_WC && IsWcPair(stb, pb)) || (a == DP_U_GU && IsGuPair(stb, pb))) {
      energy = base00 + right_unpaired;
      if (energy <= delta) exps.push_back({.delta = energy, .idx0 = t04::DpIndex(st, piv, DP_P)});

      energy = base00 + right_paired;
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0 = t04::DpIndex(st, piv, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U)});
    }

    // The rest of the cases are for U and U2.
    if (a != DP_U && a != DP_U2) continue;

    if (em_->cfg.ctd == erg::EnergyCfg::Ctd::ALL || em_->cfg.ctd == erg::EnergyCfg::Ctd::NO_COAX) {
      // (   )3<   > 3' - U, U2
      energy = base01 + em_->dangle3[pl1b][pb][stb] + em_->PfUnpaired(piv);
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired,
            .idx0 = t04::DpIndex(st, piv - 1, DP_P),
            .ctd0{st, CTD_3_DANGLE}});

      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0 = t04::DpIndex(st, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U),
            .ctd0{st, CTD_3_DANGLE}});

      // 5(   )<   > 5' - U, U2
      energy = base10 + em_->dangle5[pb][stb][st1b] + em_->PfUnpaired(st);
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired,
            .idx0 = t04::DpIndex(st + 1, piv, DP_P),
            .ctd0{st + 1, CTD_5_DANGLE}});

      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0 = t04::DpIndex(st + 1, piv, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U),
            .ctd0{st + 1, CTD_5_DANGLE}});

      // .(   ).<   > Terminal mismatch - U, U2
      energy =
          base11 + em_->terminal[pl1b][pb][stb][st1b] + em_->PfUnpaired(st) + em_->PfUnpaired(piv);
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired,
            .idx0 = t04::DpIndex(st + 1, piv - 1, DP_P),
            .ctd0{st + 1, CTD_MISMATCH}});

      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0 = t04::DpIndex(st + 1, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U),
            .ctd0{st + 1, CTD_MISMATCH}});
    }

    if (em_->cfg.ctd == erg::EnergyCfg::Ctd::ALL) {
      // .(   ).<(   ) > Left coax - U, U2
      energy = base11 + em_->MismatchCoaxial(pl1b, pb, stb, st1b) + em_->PfUnpaired(st) +
          em_->PfUnpaired(piv);
      if (energy + dp[piv + 1][en][DP_U_WC] <= delta)
        exps.push_back({
            .delta = energy + dp[piv + 1][en][DP_U_WC],
            .idx0 = t04::DpIndex(st + 1, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U_WC),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{piv + 1, CTD_LCOAX_WITH_PREV},
        });
      if (energy + dp[piv + 1][en][DP_U_GU] <= delta)
        exps.push_back({
            .delta = energy + dp[piv + 1][en][DP_U_GU],
            .idx0 = t04::DpIndex(st + 1, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U_GU),
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{piv + 1, CTD_LCOAX_WITH_PREV},
        });

      // (   )<.(   ). > Right coax forward - U, U2
      energy = base00 + dp[piv + 1][en][DP_U_RC];
      if (energy <= delta)
        exps.push_back({
            .delta = energy,
            .idx0 = t04::DpIndex(st, piv, DP_P),
            .idx1 = t04::DpIndex(piv + 1, en, DP_U_RC),
            .ctd0{st, CTD_RC_WITH_NEXT},
            .ctd1{piv + 2, CTD_RC_WITH_PREV},
        });

      // (   )(<   ) > Flush coax - U, U2
      energy = base01 + em_->stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC];
      if (energy <= delta)
        exps.push_back({
            .delta = energy,
            .idx0 = t04::DpIndex(st, piv - 1, DP_P),
            .idx1 = t04::DpIndex(piv, en, DP_U_WC),
            .ctd0{st, CTD_FCOAX_WITH_NEXT},
            .ctd1{piv, CTD_FCOAX_WITH_PREV},
        });

      if ((IsGu(pb))) {
        energy = base01 + em_->stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU];
        if (energy <= delta)
          exps.push_back({
              .delta = energy,
              .idx0 = t04::DpIndex(st, piv - 1, DP_P),
              .idx1 = t04::DpIndex(piv, en, DP_U_GU),
              .ctd0{st, CTD_FCOAX_WITH_NEXT},
              .ctd1{piv, CTD_FCOAX_WITH_PREV},
          });
      }
    }
  }

  return exps;
}

std::vector<Expansion> SuboptPersistent::PenultimateExpansions(
    int st, int en, int length, Energy delta) const {
  const auto& nostack = dp_.nostack;
  const auto& penult = dp_.penult;
  std::vector<Expansion> exps;
  Energy energy = ZERO_E;

  const auto bulge_left = em_->Bulge(r_, st, en, st + 2, en - 1);
  const auto bulge_right = em_->Bulge(r_, st, en, st + 1, en - 2);

  auto none = em_->stack[r_[st]][r_[st + 1]][r_[en - 1]][r_[en]] + em_->PfPaired(st, en) -
      penult[st][en][length];
  if (length == 2) {
    energy = none + nostack[st + 1][en - 1] +
        em_->penultimate_stack[r_[st]][r_[st + 1]][r_[en - 1]][r_[en]];
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
        em_->penultimate_stack[r_[st]][r_[st + 2]][r_[en - 1]][r_[en]];
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
        em_->penultimate_stack[r_[st]][r_[st + 1]][r_[en - 2]][r_[en]];
    if (energy <= delta)
      exps.push_back({.delta = energy, .idx0 = NoStackIndex(st + 1, en - 2), .pair{st, en}});
  }

  energy = right + penult[st + 1][en - 2][length - 1];
  if (energy <= delta)
    exps.push_back(
        {.delta = energy, .idx0 = PenultimateIndex(st + 1, en - 2, length - 1), .pair{st, en}});

  return exps;
}

}  // namespace mrna::md::t22
