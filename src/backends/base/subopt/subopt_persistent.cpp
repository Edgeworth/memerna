// Copyright 2023 Eliot Courtney.
#include "backends/base/subopt/subopt_persistent.h"

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
#include "model/energy.h"
#include "model/secondary.h"
#include "util/error.h"

namespace mrna::md::base {

namespace {

// Check the time every `CHECK_TIME_FREQ` iterations.
constexpr int CHECK_TIME_FREQ = 10000;

}  // namespace

template <bool UseLru>
SuboptPersistent<UseLru>::SuboptPersistent(Primary r, Model::Ptr m, DpState dp, SuboptCfg cfg)
    : r_(std::move(r)), m_(std::move(m)), pc_(Primary(r_), m_), dp_(std::move(dp)), cfg_(cfg),
      cache_(r_, DpIndex::MaxLinearIndex(r_.size())) {}

template <bool UseLru>
int SuboptPersistent<UseLru>::Run(const SuboptCallback& fn) {
  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false, true},
      .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::D2,
          erg::EnergyCfg::Ctd::NONE},
  };
  support.VerifySupported(funcname(), m_->cfg());
  m_->pf.Verify(r_);

  spdlog::debug("base {} with cfg {}", funcname(), m_->cfg());

  q_.clear();
  pq_ = {};  // priority queue has no clear method
  q_.reserve(r_.size() * r_.size());

  const DpIndex start_idx{0, -1, EXT};
  const Energy mfe = dp_.Index(start_idx);
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

    // TODO(2): We could expose an API that just has the index to avoid the cost of materialising
    // the full structure here.
    num_strucs++;
    GenerateResult(idx);
    res_.energy = delta + mfe;
    fn(res_);
  }

  return num_strucs;
}

template <bool UseLru>
std::pair<Energy, int> SuboptPersistent<UseLru>::RunInternal() {
  while (!pq_.empty()) {
    auto [neg_delta, idx] = pq_.top();
    pq_.pop();
    auto& s = q_[idx];

    if (s.to_expand.st == -1) {
      // At a terminal state - this is a bit different to the other subopt implementations because
      // nodes with no expansions are put onto the queue to make sure energy is properly ordered.
      return {-neg_delta, idx};
    }

    const auto& exps = GetExpansion(s.to_expand);
    assert(!exps.empty());  // Must produce at least one expansion: {-1, -1, -1}.

    const auto& exp = exps[s.expand_idx];
    const bool has_unexpanded = exp.idx1.st != -1;
    // If we have an unexpanded DpIndex, tell the child to look at us first. Otherwise, progress
    // to the next one.
    Node ns = {.parent_idx = idx,
        .parent_expand_idx = s.expand_idx,
        .unexpanded_idx = has_unexpanded ? idx : s.unexpanded_idx,
        .unexpanded_expand_idx = has_unexpanded ? s.expand_idx : s.unexpanded_expand_idx,
        .expand_idx = 0,
        .to_expand = exp.idx0};
    // Update this node for the next expansion.
    s.expand_idx++;

    // If we still have expansions to process, update the energy of the current node with what
    // the best we could do is with the next (worse) expansion.
    if (s.expand_idx != static_cast<int>(exps.size()))
      pq_.emplace(neg_delta + exp.delta - exps[s.expand_idx].delta, idx);

    if (exp.idx0.st == -1 && s.unexpanded_idx != -1) {
      // Ran out of expansions at this node but we still have unexpanded DpIndexes to process.
      assert(!has_unexpanded);
      assert(!exp.ctd0.IsValid() && !exp.ctd1.IsValid());

      // Grab the next unexpanded DpIndex.
      const auto& unexpanded_exp =
          GetExpansion(q_[s.unexpanded_idx].to_expand)[s.unexpanded_expand_idx];
      assert(unexpanded_exp.idx1.st != -1);  // There should be a valid unexpanded DpIndex here.

      // Set the next unexpanded index to the next one in the path up the tree.
      ns.unexpanded_expand_idx = q_[ns.unexpanded_idx].unexpanded_expand_idx;
      ns.unexpanded_idx = q_[ns.unexpanded_idx].unexpanded_idx;
      ns.to_expand = unexpanded_exp.idx1;
    }

    // We want to explore nodes with the same energy depth-first, so we reach a complete structure
    // as fast as possible. Doing this depth first guarantees linear memory usage in the number of
    // structures produced. We don't need to explicilty store a depth because the increasing index
    // functions as a depth counter. If we stored a depth value, the only advantage we would get is
    // starting from a lower internal node when we 'switch' to a new place in the tree to generate a
    // new structure, but the extra memory usage is not worth it.
    pq_.emplace(neg_delta, static_cast<int>(q_.size()));
    // This is the only modification to `q_`, so access to `s` is valid until here.
    q_.push_back(ns);
  }
  return {ZERO_E, -1};
}

template <bool UseLru>
void SuboptPersistent<UseLru>::GenerateResult(int idx) {
  res_.tb.s.reset(r_.size());
  res_.tb.ctd.reset(r_.size());

  int expand_idx = q_[idx].parent_expand_idx;
  idx = q_[idx].parent_idx;
  while (idx != -1) {
    const auto to_expand = q_[idx].to_expand;
    // Get the expansion that generated the previous node.
    const auto& exp = GetExpansion(to_expand)[expand_idx];

    if (to_expand.en != -1 && to_expand.a == DP_P) {
      res_.tb.s[to_expand.st] = to_expand.en;
      res_.tb.s[to_expand.en] = to_expand.st;
    }
    exp.ctd0.MaybeApply(res_.tb.ctd);
    exp.ctd1.MaybeApply(res_.tb.ctd);

    expand_idx = q_[idx].parent_expand_idx;
    idx = q_[idx].parent_idx;
  }
}

template <bool UseLru>
std::vector<Expansion> SuboptPersistent<UseLru>::GenerateExpansions(
    const DpIndex& to_expand, Energy delta) const {
  const int N = static_cast<int>(r_.size());
  const int st = to_expand.st;
  int en = to_expand.en;
  const int a = to_expand.a;
  std::vector<Expansion> exps;
  // Temporary variable to hold energy calculations.
  Energy energy = ZERO_E;
  // Exterior loop
  if (en == -1) {
    if (a == EXT) {
      // Base case: do nothing.
      if (st == N)
        exps.push_back({.delta = ZERO_E});
      else
        // Case: No pair starting here (for EXT only)
        exps.push_back({.delta = dp_.ext[st + 1][EXT] + m_->pf.Unpaired(st) - dp_.ext[st][a],
            .idx0{st + 1, -1, EXT}});
    }
    for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = r_[st];
      const auto st1b = r_[st + 1];
      const auto enb = r_[en];
      const auto en1b = r_[en - 1];
      const auto base00 = dp_.dp[st][en][DP_P] + m_->AuGuPenalty(stb, enb) - dp_.ext[st][a];
      const auto base01 = dp_.dp[st][en - 1][DP_P] + m_->AuGuPenalty(stb, en1b) - dp_.ext[st][a];
      const auto base10 = dp_.dp[st + 1][en][DP_P] + m_->AuGuPenalty(st1b, enb) - dp_.ext[st][a];
      const auto base11 =
          dp_.dp[st + 1][en - 1][DP_P] + m_->AuGuPenalty(st1b, en1b) - dp_.ext[st][a];

      // (   )<.( * ). > Right coax backward
      if (m_->cfg().UseCoaxialStacking() && a == EXT_RC) {
        energy = base11 + m_->MismatchCoaxial(en1b, enb, stb, st1b) + m_->pf.Unpaired(st) +
            m_->pf.Unpaired(en) + dp_.ext[en + 1][EXT];
        // We don't set ctds here, since we already set them in the forward case.
        if (energy <= delta)
          exps.push_back({.delta = energy, .idx0{en + 1, -1, EXT}, .idx1{st + 1, en - 1, DP_P}});
      }

      // EXT_RC is only for the above case.
      if (a == EXT_RC) continue;

      // Cases for EXT, EXT_WC, EXT_GU.
      // (   )<   >
      // If we are at EXT then this is unused.
      energy = base00 + dp_.ext[en + 1][EXT];
      Ctd val_ctd = CTD_UNUSED;

      if (m_->cfg().UseD2()) {
        // Note that D2 can overlap with anything.
        if (st != 0 && en != N - 1) {
          // (   )<   > Terminal mismatch
          energy += m_->terminal[enb][r_[en + 1]][r_[st - 1]][stb];
          val_ctd = CTD_MISMATCH;
        } else if (en != N - 1) {
          // (   )<3   > 3'
          energy += m_->dangle3[enb][r_[en + 1]][stb];
          val_ctd = CTD_3_DANGLE;
        } else if (st != 0) {
          // 5(   )<   > 5'
          energy += m_->dangle5[enb][r_[st - 1]][stb];
          val_ctd = CTD_5_DANGLE;
        }
      }

      if (energy <= delta) {
        if (a == EXT)
          exps.push_back(
              {.delta = energy, .idx0{en + 1, -1, EXT}, .idx1{st, en, DP_P}, .ctd0{st, val_ctd}});

        // (   )<   >
        // If we are at EXT_WC or EXT_GU, the CTDs for this have already have been set from a
        // coaxial stack.
        if ((a == EXT_WC && IsWcPair(stb, enb)) || (a == EXT_GU && IsGuPair(stb, enb)))
          exps.push_back({.delta = energy, .idx0{en + 1, -1, EXT}, .idx1{st, en, DP_P}});
      }

      // Everything after this is only for EXT.
      if (a != EXT) continue;

      if (m_->cfg().UseDangleMismatch()) {
        // (   )3<   > 3'
        energy = base01 + m_->dangle3[en1b][enb][stb] + m_->pf.Unpaired(en) + dp_.ext[en + 1][EXT];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{en + 1, -1, EXT},
              .idx1{st, en - 1, DP_P},
              .ctd0{st, CTD_3_DANGLE}});

        // 5(   )<   > 5'
        energy = base10 + m_->dangle5[enb][stb][st1b] + m_->pf.Unpaired(st) + dp_.ext[en + 1][EXT];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{en + 1, -1, EXT},
              .idx1{st + 1, en, DP_P},
              .ctd0{st + 1, CTD_5_DANGLE}});

        // .(   ).<   > Terminal mismatch
        energy = base11 + m_->terminal[en1b][enb][stb][st1b] + m_->pf.Unpaired(st) +
            m_->pf.Unpaired(en) + dp_.ext[en + 1][EXT];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{en + 1, -1, EXT},
              .idx1{st + 1, en - 1, DP_P},
              .ctd0{st + 1, CTD_MISMATCH}});
      }

      if (m_->cfg().UseCoaxialStacking()) {
        if (en < N - 1) {
          // .(   ).<(   ) > Left coax
          energy = base11 + m_->MismatchCoaxial(en1b, enb, stb, st1b) + m_->pf.Unpaired(st) +
              m_->pf.Unpaired(en);
          if (energy + dp_.ext[en + 1][EXT_GU] <= delta)
            exps.push_back({.delta = energy + dp_.ext[en + 1][EXT_GU],
                .idx0{en + 1, -1, EXT_GU},
                .idx1{st + 1, en - 1, DP_P},
                .ctd0{en + 1, CTD_LCOAX_WITH_PREV},
                .ctd1{st + 1, CTD_LCOAX_WITH_NEXT}});
          if (energy + dp_.ext[en + 1][EXT_WC] <= delta)
            exps.push_back({.delta = energy + dp_.ext[en + 1][EXT_WC],
                .idx0{en + 1, -1, EXT_WC},
                .idx1{st + 1, en - 1, DP_P},
                .ctd0{en + 1, CTD_LCOAX_WITH_PREV},
                .ctd1{st + 1, CTD_LCOAX_WITH_NEXT}});
        }

        if (en < N - 2) {
          // (   )<.(   ). > Right coax forward
          energy = base00 + dp_.ext[en + 1][EXT_RC];
          if (energy <= delta)
            exps.push_back({.delta = energy,
                .idx0{en + 1, -1, EXT_RC},
                .idx1{st, en, DP_P},
                .ctd0{en + 2, CTD_RC_WITH_PREV},
                .ctd1{st, CTD_RC_WITH_NEXT}});
        }

        // (   )(<   ) > Flush coax
        energy = base01 + m_->stack[en1b][enb][WcPair(enb)][stb] + dp_.ext[en][EXT_WC];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{en, -1, EXT_WC},
              .idx1{st, en - 1, DP_P},
              .ctd0{en, CTD_FCOAX_WITH_PREV},
              .ctd1{st, CTD_FCOAX_WITH_NEXT}});

        if (IsGu(enb)) {
          energy = base01 + m_->stack[en1b][enb][GuPair(enb)][stb] + dp_.ext[en][EXT_GU];
          if (energy <= delta)
            exps.push_back({.delta = energy,
                .idx0{en, -1, EXT_GU},
                .idx1{st, en - 1, DP_P},
                .ctd0{en, CTD_FCOAX_WITH_PREV},
                .ctd1{st, CTD_FCOAX_WITH_NEXT}});
        }
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
        energy = pc_.TwoLoop(st, en, ist, ien) + dp_.dp[ist][ien][DP_P] - dp_.dp[st][en][a];
        if (energy <= delta) exps.push_back({.delta = energy, .idx0{ist, ien, DP_P}});
      }
    }

    // Hairpin loop
    energy = pc_.Hairpin(st, en) - dp_.dp[st][en][a];
    if (energy <= delta) exps.push_back({.delta = energy});

    auto base_and_branch =
        pc_.augubranch[stb][enb] + m_->pf.Paired(st, en) + m_->multiloop_a - dp_.dp[st][en][a];

    // (<   ><    >)
    energy = base_and_branch + dp_.dp[st + 1][en - 1][DP_U2];
    Ctd val_ctd = CTD_UNUSED;
    if (m_->cfg().UseD2()) {
      // D2 can overlap terminal mismatches with anything.
      // (<   ><   >) Terminal mismatch
      energy += m_->terminal[stb][st1b][en1b][enb];
      val_ctd = CTD_MISMATCH;
    }
    if (energy <= delta)
      exps.push_back({.delta = energy, .idx0{st + 1, en - 1, DP_U2}, .ctd0{en, val_ctd}});

    if (m_->cfg().UseDangleMismatch()) {
      // (3<   ><   >) 3'
      energy = base_and_branch + dp_.dp[st + 2][en - 1][DP_U2] + m_->dangle3[stb][st1b][enb] +
          m_->pf.Unpaired(st + 1) + m_->multiloop_c;
      if (energy <= delta)
        exps.push_back({.delta = energy, .idx0{st + 2, en - 1, DP_U2}, .ctd0{en, CTD_3_DANGLE}});
      // (<   ><   >5) 5'
      energy = base_and_branch + dp_.dp[st + 1][en - 2][DP_U2] + m_->dangle5[stb][en1b][enb] +
          m_->pf.Unpaired(en - 1) + m_->multiloop_c;
      if (energy <= delta)
        exps.push_back({.delta = energy, .idx0{st + 1, en - 2, DP_U2}, .ctd0{en, CTD_5_DANGLE}});
      // (.<   ><   >.) Terminal mismatch
      energy = base_and_branch + dp_.dp[st + 2][en - 2][DP_U2] +
          m_->terminal[stb][st1b][en1b][enb] + m_->pf.Unpaired(st + 1) + m_->pf.Unpaired(en - 1) +
          2 * m_->multiloop_c;
      if (energy <= delta)
        exps.push_back({.delta = energy, .idx0{st + 2, en - 2, DP_U2}, .ctd0{en, CTD_MISMATCH}});
    }

    if (m_->cfg().UseCoaxialStacking()) {
      const auto outer_coax = m_->MismatchCoaxial(stb, st1b, en1b, enb) + m_->pf.Unpaired(st + 1) +
          m_->pf.Unpaired(en - 1) + 2 * m_->multiloop_c;
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
            m_->pf.Unpaired(st + 1) + m_->pf.Unpaired(piv) + 2 * m_->multiloop_c;
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st + 2, piv - 1, DP_P},
              .idx1{piv + 1, en - 1, DP_U},
              .ctd0{st + 2, CTD_RC_WITH_PREV},
              .ctd1{en, CTD_RC_WITH_NEXT}});

        // (   .(   ).) Right inner coax
        energy = base_and_branch + dp_.dp[st + 1][piv][DP_U] + pc_.augubranch[pr1b][en2b] +
            dp_.dp[piv + 2][en - 2][DP_P] + m_->MismatchCoaxial(en2b, en1b, prb, pr1b) +
            m_->pf.Unpaired(piv + 1) + m_->pf.Unpaired(en - 1) + 2 * m_->multiloop_c;
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

  // Left unpaired. Either DP_U or DP_U2.
  if (st + 1 < en && (a == DP_U || a == DP_U2)) {
    energy = dp_.dp[st + 1][en][a] + m_->pf.Unpaired(st) + m_->multiloop_c - dp_.dp[st][en][a];
    if (energy <= delta) exps.push_back({.delta = energy, .idx0{st + 1, en, a}});
  }

  // Pair here.
  for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
    //   (   .   )<   (
    // stb pl1b pb   pr1b
    auto pb = r_[piv];
    auto pl1b = r_[piv - 1];
    // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
    auto base00 = dp_.dp[st][piv][DP_P] + pc_.augubranch[stb][pb] - dp_.dp[st][en][a];
    auto base01 = dp_.dp[st][piv - 1][DP_P] + pc_.augubranch[stb][pl1b] - dp_.dp[st][en][a];
    auto base10 = dp_.dp[st + 1][piv][DP_P] + pc_.augubranch[st1b][pb] - dp_.dp[st][en][a];
    auto base11 = dp_.dp[st + 1][piv - 1][DP_P] + pc_.augubranch[st1b][pl1b] - dp_.dp[st][en][a];

    const auto right_paired = dp_.dp[piv + 1][en][DP_U];
    // This is only usable if a != DP_U2 since this leaves everything unpaired.
    const auto right_unpaired = m_->pf.UnpairedCum(piv + 1, en) + (en - piv - 1) * m_->multiloop_c;

    // Check a == U_RC:
    // (   )<.( ** ). > Right coax backward
    if (m_->cfg().UseCoaxialStacking() && a == DP_U_RC) {
      energy = base11 + m_->MismatchCoaxial(pl1b, pb, stb, st1b) + m_->pf.Unpaired(st) +
          m_->pf.Unpaired(piv) + 2 * m_->multiloop_c;
      // Our ctds will have already been set by now.
      if (energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired, .idx0{st + 1, piv - 1, DP_P}});
      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0{st + 1, piv - 1, DP_P},
            .idx1{piv + 1, en, DP_U}});
    }

    // DP_U_RC is only the above case.
    if (a == DP_U_RC) continue;

    // From here on, a must be U, U2, U_WC, or U_GU.

    // (   )<   > - U, U2, U_WC?, U_GU?
    energy = base00;
    Ctd val_ctd = CTD_UNUSED;
    if (m_->cfg().UseD2()) {
      // Note that D2 can overlap with anything.
      if (st != 0 && piv != N - 1) {
        // (   )<   > Terminal mismatch - U, U2
        energy += m_->terminal[pb][r_[piv + 1]][r_[st - 1]][stb];
        val_ctd = CTD_MISMATCH;
      } else if (piv != N - 1) {
        // (   )<3   > 3' - U, U2
        energy += m_->dangle3[pb][r_[piv + 1]][stb];
        val_ctd = CTD_3_DANGLE;
      } else if (st != 0) {
        // 5(   )<   > 5' - U, U2
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

    if ((a == DP_U_WC && IsWcPair(stb, pb)) || (a == DP_U_GU && IsGuPair(stb, pb))) {
      if (energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired, .idx0{st, piv, DP_P}});
      if (energy + right_paired <= delta)
        exps.push_back(
            {.delta = energy + right_paired, .idx0{st, piv, DP_P}, .idx1{piv + 1, en, DP_U}});
    }

    // The rest of the cases are for U and U2.
    if (a != DP_U && a != DP_U2) continue;

    if (m_->cfg().UseDangleMismatch()) {
      // (   )3<   > 3' - U, U2
      energy = base01 + m_->dangle3[pl1b][pb][stb] + m_->pf.Unpaired(piv) + m_->multiloop_c;
      // Can only let the rest be unpaired if we only need one branch, i.e. DP_U not DP_U2.
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back(
            {.delta = energy + right_unpaired, .idx0{st, piv - 1, DP_P}, .ctd0{st, CTD_3_DANGLE}});
      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0{st, piv - 1, DP_P},
            .idx1{piv + 1, en, DP_U},
            .ctd0{st, CTD_3_DANGLE}});

      // 5(   )<   > 5' - U, U2
      energy = base10 + m_->dangle5[pb][stb][st1b] + m_->pf.Unpaired(st) + m_->multiloop_c;
      if (a == DP_U && energy + right_unpaired <= delta)
        exps.push_back({.delta = energy + right_unpaired,
            .idx0{st + 1, piv, DP_P},
            .ctd0{st + 1, CTD_5_DANGLE}});
      if (energy + right_paired <= delta)
        exps.push_back({.delta = energy + right_paired,
            .idx0{st + 1, piv, DP_P},
            .idx1{piv + 1, en, DP_U},
            .ctd0{st + 1, CTD_5_DANGLE}});

      // .(   ).<   > Terminal mismatch - U, U2
      energy = base11 + m_->terminal[pl1b][pb][stb][st1b] + m_->pf.Unpaired(st) +
          m_->pf.Unpaired(piv) + 2 * m_->multiloop_c;
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

    if (m_->cfg().UseCoaxialStacking()) {
      // .(   ).<(   ) > Left coax - U, U2
      energy = base11 + m_->MismatchCoaxial(pl1b, pb, stb, st1b) + m_->pf.Unpaired(st) +
          m_->pf.Unpaired(piv) + 2 * m_->multiloop_c;
      if (energy + dp_.dp[piv + 1][en][DP_U_WC] <= delta)
        exps.push_back({.delta = energy + dp_.dp[piv + 1][en][DP_U_WC],
            .idx0{st + 1, piv - 1, DP_P},
            .idx1{piv + 1, en, DP_U_WC},
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{piv + 1, CTD_LCOAX_WITH_PREV}});
      if (energy + dp_.dp[piv + 1][en][DP_U_GU] <= delta)
        exps.push_back({.delta = energy + dp_.dp[piv + 1][en][DP_U_GU],
            .idx0{st + 1, piv - 1, DP_P},
            .idx1{piv + 1, en, DP_U_GU},
            .ctd0{st + 1, CTD_LCOAX_WITH_NEXT},
            .ctd1{piv + 1, CTD_LCOAX_WITH_PREV}});

      // (   )<.(   ). > Right coax forward - U, U2
      energy = base00 + dp_.dp[piv + 1][en][DP_U_RC];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0{st, piv, DP_P},
            .idx1{piv + 1, en, DP_U_RC},
            .ctd0{st, CTD_RC_WITH_NEXT},
            .ctd1{piv + 2, CTD_RC_WITH_PREV}});

      // (   )(<   ) > Flush coax - U, U2
      energy = base01 + m_->stack[pl1b][pb][WcPair(pb)][stb] + dp_.dp[piv][en][DP_U_WC];
      if (energy <= delta)
        exps.push_back({.delta = energy,
            .idx0{st, piv - 1, DP_P},
            .idx1{piv, en, DP_U_WC},
            .ctd0{st, CTD_FCOAX_WITH_NEXT},
            .ctd1{piv, CTD_FCOAX_WITH_PREV}});

      if (IsGu(pb)) {
        energy = base01 + m_->stack[pl1b][pb][GuPair(pb)][stb] + dp_.dp[piv][en][DP_U_GU];
        if (energy <= delta)
          exps.push_back({.delta = energy,
              .idx0{st, piv - 1, DP_P},
              .idx1{piv, en, DP_U_GU},
              .ctd0{st, CTD_FCOAX_WITH_NEXT},
              .ctd1{piv, CTD_FCOAX_WITH_PREV}});
      }
    }
  }

  return exps;
}

template class SuboptPersistent<false>;
template class SuboptPersistent<true>;

}  // namespace mrna::md::base
