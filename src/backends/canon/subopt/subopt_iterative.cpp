// Copyright 2026 Eliot Courtney.
#include "backends/canon/subopt/subopt_iterative.h"

#include <fmt/core.h>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "api/trace/trace.h"
#include "backends/common/static.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "util/error.h"
#include "util/log.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT, bool UseLru>
bool SuboptIterative<CFG, PseudofreeCfgT, UseLru>::IsSupported(
    const PseudofreeCfgT& /*pf*/, const SuboptCfg& /*subopt_cfg*/, std::string* reason) {
  if constexpr (CFG.lonely_pairs != erg::EnergyCfg::LonelyPairs::HEURISTIC &&
      CFG.lonely_pairs != erg::EnergyCfg::LonelyPairs::ON) {
    if (reason) *reason = fmt::format("lonely_pairs={} not supported", CFG.lonely_pairs);
    return false;
  }
  return true;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT, bool UseLru>
SuboptIterative<CFG, PseudofreeCfgT, UseLru>::SuboptIterative(
    Primary r, Model::Ptr m, DpState dp, PseudofreeCfgT pf, SuboptCfg subopt_cfg)
    : r_(std::move(r)), m_(std::move(m)), pf_(std::move(pf)), pc_(Primary(r_), m_, pf_),
      dp_(std::move(dp)), subopt_cfg_(subopt_cfg), cache_(r_, DpIndex::MaxLinearIndex(r_.size())) {}

template <erg::EnergyCfg CFG, class PseudofreeCfgT, bool UseLru>
int64_t SuboptIterative<CFG, PseudofreeCfgT, UseLru>::Run(const SuboptCallback& fn) {
  res_ = SuboptResult(ZERO_E, trace::TraceResult(Secondary(r_.size()), Ctds(r_.size())));
  q_.reserve(r_.size());

  std::string reason;
  verify(IsSupported(pf_, subopt_cfg_, &reason), "{} does not support the given configuration: {}",
      funcname(), reason);
  pf_.Verify(r_);

  logdebug("canon {} with {}, {}, {}", funcname(), CFG, subopt_cfg_, pf_);

  if (subopt_cfg_.sorted || subopt_cfg_.strucs != SuboptCfg::MAX_STRUCTURES ||
      subopt_cfg_.time_secs >= 0.0) {
    int64_t count = 0;
    Energy delta = ZERO_E;
    auto start_time = std::chrono::steady_clock::now();
    while (count < subopt_cfg_.strucs && delta != MAX_E && delta <= subopt_cfg_.delta) {
      if (subopt_cfg_.time_secs >= 0.0) {
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
            std::chrono::steady_clock::now() - start_time);
        if (elapsed.count() >= subopt_cfg_.time_secs) break;
      }
      auto res = RunInternal(fn, delta, /*exact_energy=*/true, subopt_cfg_.strucs - count);
      count += res.count;
      if (!res.keep_going) break;
      delta = res.next_seen;
    }
    return count;
  }
  return RunInternal(fn, subopt_cfg_.delta, /*exact_energy=*/false, subopt_cfg_.strucs).count;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT, bool UseLru>
typename SuboptIterative<CFG, PseudofreeCfgT, UseLru>::RunInternalResult
SuboptIterative<CFG, PseudofreeCfgT, UseLru>::RunInternal(
    const SuboptCallback& fn, Energy delta, bool exact_energy, int64_t max) {
  int64_t count = 0;
  Energy next_seen = MAX_E;
  Energy energy = ZERO_E;
  q_.clear();
  unexpanded_.clear();
  q_.push_back({.expand_idx = 0, .to_expand{0, INVALID_INDEX, EXT}, .should_unexpand = false});
  while (!q_.empty()) {
    auto& s = q_.back();
    assert(s.to_expand.IsValid());

    const auto& exps = GetExpansion(s.to_expand);
    const auto num_exps = As<ExpansionIndex>(exps.size());
    assert(!exps.empty());

    if (s.expand_idx != 0) {
      const auto& pexp = exps[s.expand_idx - 1];
      pexp.ctd0.MaybeRemove(res_.tb.ctd);
      pexp.ctd1.MaybeRemove(res_.tb.ctd);
      pexp.ctd2.MaybeRemove(res_.tb.ctd);
      if (pexp.idx2.IsValid()) unexpanded_.pop_back();
      if (pexp.idx1.IsValid()) unexpanded_.pop_back();
      energy -= pexp.delta;
    }

    if (s.expand_idx != num_exps && exps[s.expand_idx].delta + energy > delta)
      next_seen = std::min(next_seen, exps[s.expand_idx].delta + energy);

    if (s.expand_idx == num_exps || exps[s.expand_idx].delta + energy > delta) {
      if (!s.to_expand.IsExternal() && s.to_expand.a == DP_P)
        res_.tb.s[s.to_expand.st] = res_.tb.s[s.to_expand.en] = INVALID_INDEX;
      if (s.should_unexpand) unexpanded_.push_back(s.to_expand);
      q_.pop_back();
      continue;
    }

    const auto& exp = exps[s.expand_idx++];
    Node ns = {.expand_idx = 0, .to_expand = exp.idx0, .should_unexpand = false};
    energy += exp.delta;

    if (!exp.idx0.IsValid()) {
      assert(!exp.idx1.IsValid());
      assert(!exp.idx2.IsValid());
      assert(!exp.ctd0.IsValid() && !exp.ctd1.IsValid() && !exp.ctd2.IsValid());

      if (unexpanded_.empty()) {
        if (!exact_energy || energy == delta) {
          res_.energy = energy + dp_.ext[0][EXT];
          ++count;
          if (!fn(res_)) return {.count = count, .next_seen = next_seen, .keep_going = false};
          if (count == max) return {.count = count, .next_seen = CAP_E};
        }
        continue;
      }
      ns.to_expand = unexpanded_.back();
      unexpanded_.pop_back();
      ns.should_unexpand = true;
    } else {
      assert(!exp.idx2.IsValid() || exp.idx1.IsValid());
      exp.ctd0.MaybeApply(res_.tb.ctd);
      exp.ctd1.MaybeApply(res_.tb.ctd);
      exp.ctd2.MaybeApply(res_.tb.ctd);
      if (exp.idx2.IsValid()) unexpanded_.push_back(exp.idx2);
      if (exp.idx1.IsValid()) unexpanded_.push_back(exp.idx1);
    }

    if (!ns.to_expand.IsExternal() && ns.to_expand.a == DP_P) {
      res_.tb.s[ns.to_expand.st] = ns.to_expand.en;
      res_.tb.s[ns.to_expand.en] = ns.to_expand.st;
    }

    q_.push_back(ns);
  }
  assert(unexpanded_.empty() && energy == ZERO_E && res_.tb.s == Secondary(res_.tb.s.size()) &&
      res_.tb.ctd == Ctds(res_.tb.ctd.size()));
  return {.count = count, .next_seen = next_seen};
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT, bool UseLru>
std::vector<Expansion> SuboptIterative<CFG, PseudofreeCfgT, UseLru>::GenerateExpansions(
    const DpIndex& to_expand, Energy delta) const {
  const int N = r_.size();
  const int st = to_expand.st;
  int en = to_expand.en;
  const DpArrayId a = to_expand.a;
  std::vector<Expansion> exps;
  Energy energy = ZERO_E;

  if (to_expand.IsExternal()) {
    if (a == EXT) {
      if (st == N)
        exps.push_back({.delta = ZERO_E});
      else
        exps.push_back({.delta = dp_.ext[st + 1][EXT] + pf_.Unpaired(st) - dp_.ext[st][a],
            .idx0{st + 1, INVALID_INDEX, EXT}});
    }
    for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const auto stb = r_[st];
      const auto st1b = r_[st + 1];
      const auto enb = r_[en];
      const auto en1b = r_[en - 1];
      const auto base00 = dp_.dp[st][en][DP_P] + m_->AuGuPenalty(stb, enb) - dp_.ext[st][a];
      const auto base01 = dp_.dp[st][en - 1][DP_P] + m_->AuGuPenalty(stb, en1b) - dp_.ext[st][a];
      const auto base10 = dp_.dp[st + 1][en][DP_P] + m_->AuGuPenalty(st1b, enb) - dp_.ext[st][a];
      const auto base11 =
          dp_.dp[st + 1][en - 1][DP_P] + m_->AuGuPenalty(st1b, en1b) - dp_.ext[st][a];

      energy = base00 + dp_.ext[en + 1][EXT];
      Ctd val_ctd = CTD_UNUSED;
      if constexpr (CFG.UseD2()) {
        if (st != 0 && en != N - 1) {
          energy += m_->terminal[enb][r_[en + 1]][r_[st - 1]][stb];
          val_ctd = CTD_MISMATCH;
        } else if (en != N - 1) {
          energy += m_->dangle3[enb][r_[en + 1]][stb];
          val_ctd = CTD_3_DANGLE;
        } else if (st != 0) {
          energy += m_->dangle5[enb][r_[st - 1]][stb];
          val_ctd = CTD_5_DANGLE;
        }
      }
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0{en + 1, INVALID_INDEX, EXT},
            .idx1{st, en, DP_P},
            .ctd0{st, val_ctd}});

      if constexpr (CFG.UseDangleMismatch()) {
        energy = base01 + m_->dangle3[en1b][enb][stb] + pf_.Unpaired(en) + dp_.ext[en + 1][EXT];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{en + 1, INVALID_INDEX, EXT},
              .idx1{st, en - 1, DP_P},
              .ctd0{st, CTD_3_DANGLE}});

        energy = base10 + m_->dangle5[enb][stb][st1b] + pf_.Unpaired(st) + dp_.ext[en + 1][EXT];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{en + 1, INVALID_INDEX, EXT},
              .idx1{st + 1, en, DP_P},
              .ctd0{st + 1, CTD_5_DANGLE}});

        energy = base11 + m_->terminal[en1b][enb][stb][st1b] + pf_.Unpaired(st) + pf_.Unpaired(en) +
            dp_.ext[en + 1][EXT];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{en + 1, INVALID_INDEX, EXT},
              .idx1{st + 1, en - 1, DP_P},
              .ctd0{st + 1, CTD_MISMATCH}});
      }

      if constexpr (CFG.UseCoaxialStacking()) {
        // .(   ).<(   ) > Left coax
        for (int piv = en + HAIRPIN_MIN_SZ + 2; piv < N; ++piv) {
          const Base prb = r_[en + 1];
          const Base pb = r_[piv];
          energy = base11 + dp_.dp[en + 1][piv][DP_P] + m_->AuGuPenalty(prb, pb) +
              m_->MismatchCoaxial(en1b, enb, stb, st1b) + pf_.Unpaired(st) + pf_.Unpaired(en) +
              dp_.ext[piv + 1][EXT];
          if (energy <= delta)
            exps.push_back({.delta = energy,
                .idx0{piv + 1, INVALID_INDEX, EXT},
                .idx1{st + 1, en - 1, DP_P},
                .idx2{en + 1, piv, DP_P},
                .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
                .ctd1{en + 1, CTD_LCOAX_WITH_PREV}});
        }

        // (   )<.(   ). > Right coax forward
        for (int piv = en + HAIRPIN_MIN_SZ + 4; piv < N; ++piv) {
          const Base prb = r_[en + 1];
          const Base pr1b = r_[en + 2];
          const Base pl1b = r_[piv - 1];
          const Base pb = r_[piv];
          energy = base00 + dp_.dp[en + 2][piv - 1][DP_P] + m_->AuGuPenalty(pr1b, pl1b) +
              m_->MismatchCoaxial(pl1b, pb, prb, pr1b) + pf_.Unpaired(en + 1) + pf_.Unpaired(piv) +
              dp_.ext[piv + 1][EXT];
          if (energy <= delta)
            exps.push_back({.delta = energy,
                .idx0{piv + 1, INVALID_INDEX, EXT},
                .idx1{st, en, DP_P},
                .idx2{en + 2, piv - 1, DP_P},
                .ctd0{st, CTD_RC_WITH_NEXT},
                .ctd1{en + 2, CTD_RC_WITH_PREV}});
        }

        // (   )(<   ) > Flush coax
        for (int piv = en + HAIRPIN_MIN_SZ + 1; piv < N; ++piv) {
          const Base pb = r_[piv];
          energy = base01 + dp_.dp[en][piv][DP_P] + m_->AuGuPenalty(enb, pb) +
              m_->stack[en1b][enb][pb][stb] + dp_.ext[piv + 1][EXT];
          if (energy <= delta)
            exps.push_back({.delta = energy,
                .idx0{piv + 1, INVALID_INDEX, EXT},
                .idx1{st, en - 1, DP_P},
                .idx2{en, piv, DP_P},
                .ctd0{st, CTD_FCOAX_WITH_NEXT},
                .ctd1{en, CTD_FCOAX_WITH_PREV}});
        }
      }
    }
    return exps;
  }

  const auto stb = r_[st];
  const auto st1b = r_[st + 1];
  const auto st2b = r_[st + 2];
  const auto enb = r_[en];
  const auto en1b = r_[en - 1];
  const auto en2b = r_[en - 2];

  if (a == DP_P) {
    const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
    for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
      for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
        energy = pc_.TwoLoop(st, en, ist, ien) + dp_.dp[ist][ien][DP_P] - dp_.dp[st][en][a];
        if (energy <= delta) exps.push_back({.delta = energy, .idx0{ist, ien, DP_P}});
      }
    }

    energy = pc_.Hairpin(st, en) - dp_.dp[st][en][a];
    if (energy <= delta) exps.push_back({.delta = energy});

    auto base_and_branch =
        pc_.augubranch[stb][enb] + pf_.Paired(st, en) + m_->multiloop_a - dp_.dp[st][en][a];

    energy = base_and_branch + dp_.dp[st + 1][en - 1][DP_U2];
    Ctd val_ctd = CTD_UNUSED;
    if constexpr (CFG.UseD2()) {
      energy += m_->terminal[stb][st1b][en1b][enb];
      val_ctd = CTD_MISMATCH;
    }
    if (energy <= delta)
      exps.push_back({.delta = energy, .idx0{st + 1, en - 1, DP_U2}, .ctd0{en, val_ctd}});

    if constexpr (CFG.UseDangleMismatch()) {
      energy = base_and_branch + dp_.dp[st + 2][en - 1][DP_U2] + m_->dangle3[stb][st1b][enb] +
          pf_.Unpaired(st + 1) + m_->multiloop_c;
      if (energy <= delta)
        exps.push_back({.delta = energy, .idx0{st + 2, en - 1, DP_U2}, .ctd0{en, CTD_3_DANGLE}});

      energy = base_and_branch + dp_.dp[st + 1][en - 2][DP_U2] + m_->dangle5[stb][en1b][enb] +
          pf_.Unpaired(en - 1) + m_->multiloop_c;
      if (energy <= delta)
        exps.push_back({.delta = energy, .idx0{st + 1, en - 2, DP_U2}, .ctd0{en, CTD_5_DANGLE}});

      energy = base_and_branch + dp_.dp[st + 2][en - 2][DP_U2] +
          m_->terminal[stb][st1b][en1b][enb] + pf_.Unpaired(st + 1) + pf_.Unpaired(en - 1) +
          2 * m_->multiloop_c;
      if (energy <= delta)
        exps.push_back({.delta = energy, .idx0{st + 2, en - 2, DP_U2}, .ctd0{en, CTD_MISMATCH}});
    }

    if constexpr (CFG.UseCoaxialStacking()) {
      const auto outer_coax = m_->MismatchCoaxial(stb, st1b, en1b, enb) + pf_.Unpaired(st + 1) +
          pf_.Unpaired(en - 1) + 2 * m_->multiloop_c;
      for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
        const Base pl1b = r_[piv - 1];
        const Base plb = r_[piv];
        const Base prb = r_[piv + 1];
        const Base pr1b = r_[piv + 2];

        // (.(   )   .) Left outer coax - P
        energy = base_and_branch + dp_.dp[st + 2][piv][DP_P] + pc_.augubranch[st2b][plb] +
            dp_.dp[piv + 1][en - 2][DP_U] + outer_coax;
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st + 2, piv, DP_P},
              .idx1{piv + 1, en - 2, DP_U},
              .ctd0{st + 2, CTD_LCOAX_WITH_PREV},
              .ctd1{en, CTD_LCOAX_WITH_NEXT}});

        // (.   (   ).) Right outer coax
        energy = base_and_branch + dp_.dp[st + 2][piv][DP_U] + pc_.augubranch[prb][en2b] +
            dp_.dp[piv + 1][en - 2][DP_P] + outer_coax;
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st + 2, piv, DP_U},
              .idx1{piv + 1, en - 2, DP_P},
              .ctd0{piv + 1, CTD_RC_WITH_NEXT},
              .ctd1{en, CTD_RC_WITH_PREV}});

        // (.(   ).   ) Left inner coax
        energy = base_and_branch + dp_.dp[st + 2][piv - 1][DP_P] + pc_.augubranch[st2b][pl1b] +
            dp_.dp[piv + 1][en - 1][DP_U] + m_->MismatchCoaxial(pl1b, plb, st1b, st2b) +
            pf_.Unpaired(st + 1) + pf_.Unpaired(piv) + 2 * m_->multiloop_c;
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st + 2, piv - 1, DP_P},
              .idx1{piv + 1, en - 1, DP_U},
              .ctd0{st + 2, CTD_RC_WITH_PREV},
              .ctd1{en, CTD_RC_WITH_NEXT}});

        // (   .(   ).) Right inner coax
        energy = base_and_branch + dp_.dp[st + 1][piv][DP_U] + pc_.augubranch[pr1b][en2b] +
            dp_.dp[piv + 2][en - 2][DP_P] + m_->MismatchCoaxial(en2b, en1b, prb, pr1b) +
            pf_.Unpaired(piv + 1) + pf_.Unpaired(en - 1) + 2 * m_->multiloop_c;
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st + 1, piv, DP_U},
              .idx1{piv + 2, en - 2, DP_P},
              .ctd0{piv + 2, CTD_LCOAX_WITH_NEXT},
              .ctd1{en, CTD_LCOAX_WITH_PREV}});

        // ((   )   ) Left flush coax
        energy = base_and_branch + dp_.dp[st + 1][piv][DP_P] + pc_.augubranch[st1b][plb] +
            dp_.dp[piv + 1][en - 1][DP_U] + m_->stack[stb][st1b][plb][enb];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st + 1, piv, DP_P},
              .idx1{piv + 1, en - 1, DP_U},
              .ctd0{st + 1, CTD_FCOAX_WITH_PREV},
              .ctd1{en, CTD_FCOAX_WITH_NEXT}});

        // (   (   )) Right flush coax
        energy = base_and_branch + dp_.dp[st + 1][piv][DP_U] + pc_.augubranch[prb][en1b] +
            dp_.dp[piv + 1][en - 1][DP_P] + m_->stack[stb][prb][en1b][enb];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st + 1, piv, DP_U},
              .idx1{piv + 1, en - 1, DP_P},
              .ctd0{piv + 1, CTD_FCOAX_WITH_NEXT},
              .ctd1{en, CTD_FCOAX_WITH_PREV}});
      }
    }
    return exps;
  }

  if (a == DP_CX) {
    if constexpr (CFG.UseCoaxialStacking()) {
      // (   )(<   ) Flush coax - CX
      for (int piv = st + HAIRPIN_MIN_SZ + 2; piv <= en - HAIRPIN_MIN_SZ - 1; ++piv) {
        const Base pl1b = r_[piv - 1];
        const Base pb = r_[piv];
        energy = dp_.dp[st][piv - 1][DP_P] + pc_.augubranch[stb][pl1b] + dp_.dp[piv][en][DP_P] +
            pc_.augubranch[pb][enb] + m_->stack[pl1b][pb][enb][stb] - dp_.dp[st][en][a];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st, piv - 1, DP_P},
              .idx1{piv, en, DP_P},
              .ctd0{st, CTD_FCOAX_WITH_NEXT},
              .ctd1{piv, CTD_FCOAX_WITH_PREV}});
      }

      // .(   ).<(   ) Left coax - CX
      for (int piv = st + HAIRPIN_MIN_SZ + 3; piv <= en - HAIRPIN_MIN_SZ - 2; ++piv) {
        const Base pl1b = r_[piv - 1];
        const Base pb = r_[piv];
        const Base prb = r_[piv + 1];
        energy = dp_.dp[st + 1][piv - 1][DP_P] + pc_.augubranch[st1b][pl1b] +
            dp_.dp[piv + 1][en][DP_P] + pc_.augubranch[prb][enb] +
            m_->MismatchCoaxial(pl1b, pb, stb, st1b) + pf_.Unpaired(st) + pf_.Unpaired(piv) +
            2 * m_->multiloop_c - dp_.dp[st][en][a];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st + 1, piv - 1, DP_P},
              .idx1{piv + 1, en, DP_P},
              .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
              .ctd1{piv + 1, CTD_LCOAX_WITH_PREV}});
      }

      // (   )<.(   ). Right coax - CX
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en - HAIRPIN_MIN_SZ - 4; ++piv) {
        const Base pb = r_[piv];
        const Base prb = r_[piv + 1];
        const Base pr1b = r_[piv + 2];
        energy = dp_.dp[st][piv][DP_P] + pc_.augubranch[stb][pb] + dp_.dp[piv + 2][en - 1][DP_P] +
            pc_.augubranch[pr1b][en1b] + m_->MismatchCoaxial(en1b, enb, prb, pr1b) +
            pf_.Unpaired(piv + 1) + pf_.Unpaired(en) + 2 * m_->multiloop_c - dp_.dp[st][en][a];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st, piv, DP_P},
              .idx1{piv + 2, en - 1, DP_P},
              .ctd0{st, CTD_RC_WITH_NEXT},
              .ctd1{piv + 2, CTD_RC_WITH_PREV}});
      }
    }
    return exps;
  }

  if (st + 1 < en && (a == DP_U || a == DP_U2)) {
    energy = dp_.dp[st + 1][en][a] + pf_.Unpaired(st) + m_->multiloop_c - dp_.dp[st][en][a];
    if (energy <= delta) exps.push_back({.delta = energy, .idx0{st + 1, en, a}});
  }

  for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
    const auto pb = r_[piv];
    const auto pl1b = r_[piv - 1];
    const auto base00 = dp_.dp[st][piv][DP_P] + pc_.augubranch[stb][pb] - dp_.dp[st][en][a];
    const auto base01 = dp_.dp[st][piv - 1][DP_P] + pc_.augubranch[stb][pl1b] - dp_.dp[st][en][a];
    const auto base10 = dp_.dp[st + 1][piv][DP_P] + pc_.augubranch[st1b][pb] - dp_.dp[st][en][a];
    const auto base11 =
        dp_.dp[st + 1][piv - 1][DP_P] + pc_.augubranch[st1b][pl1b] - dp_.dp[st][en][a];

    const auto right_paired = dp_.dp[piv + 1][en][DP_U];
    const auto right_unpaired = pf_.UnpairedSum(piv + 1, en) + (en - piv) * m_->multiloop_c;

    energy = base00;
    Ctd val_ctd = CTD_UNUSED;
    if constexpr (CFG.UseD2()) {
      if (st != 0 && piv != N - 1) {
        energy += m_->terminal[pb][r_[piv + 1]][r_[st - 1]][stb];
        val_ctd = CTD_MISMATCH;
      } else if (piv != N - 1) {
        energy += m_->dangle3[pb][r_[piv + 1]][stb];
        val_ctd = CTD_3_DANGLE;
      } else if (st != 0) {
        energy += m_->dangle5[pb][r_[st - 1]][stb];
        val_ctd = CTD_5_DANGLE;
      }
    }
    if (a == DP_U) {
      if (energy + right_unpaired <= delta)
        exps.push_back(
            {.delta = energy + right_unpaired, .idx0{st, piv, DP_P}, .ctd0{st, val_ctd}});
      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0{st, piv, DP_P},
            .idx1{piv + 1, en, DP_U},
            .ctd0{st, val_ctd}});
    }
    if (a == DP_U2 && energy + right_paired <= delta)
      exps.push_back({.delta = energy + right_paired,
          .idx0{st, piv, DP_P},
          .idx1{piv + 1, en, DP_U},
          .ctd0{st, val_ctd}});

    if constexpr (CFG.UseDangleMismatch()) {
      energy = base01 + m_->dangle3[pl1b][pb][stb] + pf_.Unpaired(piv) + m_->multiloop_c;
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back(
            {.delta = energy + right_unpaired, .idx0{st, piv - 1, DP_P}, .ctd0{st, CTD_3_DANGLE}});
      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0{st, piv - 1, DP_P},
            .idx1{piv + 1, en, DP_U},
            .ctd0{st, CTD_3_DANGLE}});

      energy = base10 + m_->dangle5[pb][stb][st1b] + pf_.Unpaired(st) + m_->multiloop_c;
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired,
            .idx0{st + 1, piv, DP_P},
            .ctd0{st + 1, CTD_5_DANGLE}});
      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0{st + 1, piv, DP_P},
            .idx1{piv + 1, en, DP_U},
            .ctd0{st + 1, CTD_5_DANGLE}});

      energy = base11 + m_->terminal[pl1b][pb][stb][st1b] + pf_.Unpaired(st) + pf_.Unpaired(piv) +
          2 * m_->multiloop_c;
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired,
            .idx0{st + 1, piv - 1, DP_P},
            .ctd0{st + 1, CTD_MISMATCH}});
      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0{st + 1, piv - 1, DP_P},
            .idx1{piv + 1, en, DP_U},
            .ctd0{st + 1, CTD_MISMATCH}});
    }
  }

  if constexpr (CFG.UseCoaxialStacking()) {
    if (a == DP_U || a == DP_U2) {
      for (int piv = st + 2 * HAIRPIN_MIN_SZ + 3; piv <= en; ++piv) {
        const auto cx_base = dp_.dp[st][piv][DP_CX] - dp_.dp[st][en][a];
        const auto right_paired = dp_.dp[piv + 1][en][DP_U];
        const auto right_unpaired = pf_.UnpairedSum(piv + 1, en) + (en - piv) * m_->multiloop_c;
        if (cx_base + right_unpaired <= delta)
          exps.push_back({.delta = cx_base + right_unpaired, .idx0{st, piv, DP_CX}});
        if (cx_base + right_paired <= delta)
          exps.push_back(
              {.delta = cx_base + right_paired, .idx0{st, piv, DP_CX}, .idx1{piv + 1, en, DP_U}});
      }
    }
  }

  return exps;
}

#define MRNA_INSTANTIATE(CFG, PF)                 \
  template class SuboptIterative<CFG, PF, false>; \
  template class SuboptIterative<CFG, PF, true>;
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::canon
