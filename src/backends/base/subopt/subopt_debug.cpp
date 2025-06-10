// Copyright 2016 Eliot Courtney.
#include "backends/base/subopt/subopt_debug.h"

#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <utility>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "util/error.h"

namespace mrna::md::base {

SuboptDebug::SuboptDebug(Primary r, Model::Ptr m, DpState dp, SuboptCfg cfg)
    : r_(std::move(r)), m_(std::move(m)), dp_(std::move(dp)), cfg_(cfg) {}

int SuboptDebug::Run(const SuboptCallback& fn) {
  const int N = static_cast<int>(r_.size());
  verify(N < std::numeric_limits<Index>::max(), "RNA too long for suboptimal folding");

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false, true},
      .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::D2,
          erg::EnergyCfg::Ctd::NONE},
  };
  support.VerifySupported(funcname(), m_->cfg());
  m_->pf.Verify(r_);

  spdlog::debug("base {} with cfg {}", funcname(), m_->cfg());
  auto start_time = std::chrono::steady_clock::now();

  // Basic idea of suboptimal traceback is look at all possible choices from a state, and expand
  // just one of them. Fully expanding one of them means there will be no duplicates in the tree.
  // Cull the ones not inside the window or when we have more than `max_structures`.
  // We don't have to check for expanding impossible states indirectly, since they will have MAX_E,
  // be above cfg_.delta, and be instantly culled (callers use CAP_E for no energy limit).
  q_.insert({.not_yet_expanded = {{0, -1, EXT}},
      .res = SuboptResult(dp_.ext[0][EXT], trace::TraceResult(Secondary(N), Ctds(N)))});
  Node node;
  while (!q_.empty()) {
    if (cfg_.time_secs >= 0.0) {
      auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - start_time);
      if (elapsed.count() >= cfg_.time_secs) break;
    }

    node = std::move(q_.extract(q_.begin()).value());
    // Finished state.
    if (node.not_yet_expanded.empty()) {
      PruneInsert(node, &finished_);
      continue;
    }

    // If we found a non-finished node, but `finished` is full, and the worst in `finished` is
    // as good as our current node (which is the best in `q`), then we can exit.
    if (static_cast<int>(finished_.size()) >= cfg_.strucs &&
        (--finished_.end())->res.energy <= node.res.energy)
      break;

    auto to_expand = node.not_yet_expanded.back();
    node.not_yet_expanded.pop_back();
    const int st = to_expand.st;
    int en = to_expand.en;
    const int a = to_expand.a;

    // Initialise - we only make small modifications to it.
    curnode_ = node.copy();
    // Temporary variable to hold energy calculations.
    Energy energy = ZERO_E;

    // Exterior loop
    if (en == -1) {
      // We try replacing what we do at (st, a) with a bunch of different cases, so we use this
      // energy as a base.
      const Energy base_energy = node.res.energy - dp_.ext[st][a];
      if (a == EXT) {
        // Base case: do nothing.
        if (st == N)
          Expand(base_energy);
        else
          // Case: No pair starting here (for EXT only)
          Expand(base_energy + m_->pf.Unpaired(st) + dp_.ext[st + 1][EXT], {st + 1, -1, EXT});
      }
      for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        // .   .   .   (   .   .   .   )   <   >
        //           stb  st1b   en1b  enb   rem
        const auto stb = r_[st];
        const auto st1b = r_[st + 1];
        const auto enb = r_[en];
        const auto en1b = r_[en - 1];
        const auto base00 = dp_.dp[st][en][DP_P] + m_->AuGuPenalty(stb, enb);
        const auto base01 = dp_.dp[st][en - 1][DP_P] + m_->AuGuPenalty(stb, en1b);
        const auto base10 = dp_.dp[st + 1][en][DP_P] + m_->AuGuPenalty(st1b, enb);
        const auto base11 = dp_.dp[st + 1][en - 1][DP_P] + m_->AuGuPenalty(st1b, en1b);
        curnode_ = node.copy();

        // (   )<.( * ). > Right coax backward
        if (m_->cfg().UseCoaxialStacking() && a == EXT_RC) {
          energy = base_energy + base11 + m_->MismatchCoaxial(en1b, enb, stb, st1b) +
              m_->pf.Unpaired(st) + m_->pf.Unpaired(en) + dp_.ext[en + 1][EXT];
          // We don't set ctds here, since we already set them in the forward case.
          Expand(energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P});
        }

        // EXT_RC is only for the above case.
        if (a == EXT_RC) continue;

        // Cases for EXT, EXT_WC, EXT_GU.
        // (   )<   >
        energy = base_energy + base00 + dp_.ext[en + 1][EXT];
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

        if (a == EXT) Expand(energy, {en + 1, -1, EXT}, {st, en, DP_P}, {st, val_ctd});

        // (   )<   >
        // If we are at EXT_WC or EXT_GU, the CTDs for this have already have been set from a
        // coaxial stack.
        if ((a == EXT_WC && IsWcPair(stb, enb)) || (a == EXT_GU && IsGuPair(stb, enb)))
          Expand(energy, {en + 1, -1, EXT}, {st, en, DP_P});

        // Everything after this is only for EXT.
        if (a != EXT) continue;

        if (m_->cfg().UseDangleMismatch()) {
          // (   )3<   > 3'
          energy = base_energy + base01 + m_->dangle3[en1b][enb][stb] + m_->pf.Unpaired(en) +
              dp_.ext[en + 1][EXT];
          Expand(energy, {en + 1, -1, EXT}, {st, en - 1, DP_P}, {st, CTD_3_DANGLE});

          // 5(   )<   > 5'
          energy = base_energy + base10 + m_->dangle5[enb][stb][st1b] + m_->pf.Unpaired(st) +
              dp_.ext[en + 1][EXT];
          Expand(energy, {en + 1, -1, EXT}, {st + 1, en, DP_P}, {st + 1, CTD_5_DANGLE});

          // .(   ).<   > Terminal mismatch
          energy = base_energy + base11 + m_->terminal[en1b][enb][stb][st1b] + m_->pf.Unpaired(st) +
              m_->pf.Unpaired(en) + dp_.ext[en + 1][EXT];
          Expand(energy, {en + 1, -1, EXT}, {st + 1, en - 1, DP_P}, {st + 1, CTD_MISMATCH});
        }

        if (m_->cfg().UseCoaxialStacking()) {
          // (   )(<   ) > Flush coax
          energy =
              base_energy + base01 + m_->stack[en1b][enb][WcPair(enb)][stb] + dp_.ext[en][EXT_WC];
          Expand(energy, {en, -1, EXT_WC}, {st, en - 1, DP_P}, {en, CTD_FCOAX_WITH_PREV},
              {st, CTD_FCOAX_WITH_NEXT});
          if (IsGu(enb)) {
            energy =
                base_energy + base01 + m_->stack[en1b][enb][GuPair(enb)][stb] + dp_.ext[en][EXT_GU];
            Expand(energy, {en, -1, EXT_GU}, {st, en - 1, DP_P}, {en, CTD_FCOAX_WITH_PREV},
                {st, CTD_FCOAX_WITH_NEXT});
          }

          if (en < N - 1) {
            // .(   ).<(   ) > Left coax
            energy = base_energy + base11 + m_->MismatchCoaxial(en1b, enb, stb, st1b) +
                m_->pf.Unpaired(st) + m_->pf.Unpaired(en);
            Expand(energy + dp_.ext[en + 1][EXT_GU], {en + 1, -1, EXT_GU}, {st + 1, en - 1, DP_P},
                {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT});
            Expand(energy + dp_.ext[en + 1][EXT_WC], {en + 1, -1, EXT_WC}, {st + 1, en - 1, DP_P},
                {en + 1, CTD_LCOAX_WITH_PREV}, {st + 1, CTD_LCOAX_WITH_NEXT});
          }
          if (en < N - 2) {
            // (   )<.(   ). > Right coax forward
            energy = base_energy + base00 + dp_.ext[en + 1][EXT_RC];
            Expand(energy, {en + 1, -1, EXT_RC}, {st, en, DP_P}, {en + 2, CTD_RC_WITH_PREV},
                {st, CTD_RC_WITH_NEXT});
          }
        }
      }
      // Finished exterior loop, don't do anymore.
      continue;
    }

    // Subtract the minimum energy of the contribution at this node.
    const Energy base_energy = node.res.energy - dp_.dp[st][en][a];
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
      const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
      for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
        for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
          energy = base_energy + m_->TwoLoop(r_, st, en, ist, ien) + dp_.dp[ist][ien][DP_P];
          Expand(energy, {ist, ien, DP_P});
        }
      }

      // Hairpin loop
      energy = base_energy + m_->Hairpin(r_, st, en);
      Expand(energy);

      auto base_and_branch = base_energy + m_->AuGuPenalty(stb, enb) + m_->pf.Paired(st, en) +
          m_->multiloop_a + m_->multiloop_b;

      // (<   ><    >)
      energy = base_and_branch + dp_.dp[st + 1][en - 1][DP_U2];
      Ctd val_ctd = CTD_UNUSED;
      if (m_->cfg().UseD2()) {
        // D2 can overlap terminal mismatches with anything.
        // (<   ><   >) Terminal mismatch
        energy += m_->terminal[stb][st1b][en1b][enb];
        val_ctd = CTD_MISMATCH;
      }
      Expand(energy, {st + 1, en - 1, DP_U2}, {en, val_ctd});

      if (m_->cfg().UseDangleMismatch()) {
        // (3<   ><   >) 3'
        energy = base_and_branch + dp_.dp[st + 2][en - 1][DP_U2] + m_->dangle3[stb][st1b][enb] +
            m_->pf.Unpaired(st + 1);
        Expand(energy, {st + 2, en - 1, DP_U2}, {en, CTD_3_DANGLE});
        // (<   ><   >5) 5'
        energy = base_and_branch + dp_.dp[st + 1][en - 2][DP_U2] + m_->dangle5[stb][en1b][enb] +
            m_->pf.Unpaired(en - 1);
        Expand(energy, {st + 1, en - 2, DP_U2}, {en, CTD_5_DANGLE});
        // (.<   ><   >.) Terminal mismatch
        energy = base_and_branch + dp_.dp[st + 2][en - 2][DP_U2] +
            m_->terminal[stb][st1b][en1b][enb] + m_->pf.Unpaired(st + 1) + m_->pf.Unpaired(en - 1);
        Expand(energy, {st + 2, en - 2, DP_U2}, {en, CTD_MISMATCH});
      }

      if (m_->cfg().UseCoaxialStacking()) {
        const auto outer_coax = m_->MismatchCoaxial(stb, st1b, en1b, enb) +
            m_->pf.Unpaired(st + 1) + m_->pf.Unpaired(en - 1);
        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          const Base pl1b = r_[piv - 1];
          const Base plb = r_[piv];
          const Base prb = r_[piv + 1];
          const Base pr1b = r_[piv + 2];

          // (.(   )   .) Left outer coax - P
          energy = base_and_branch + dp_.dp[st + 2][piv][DP_P] + m_->multiloop_b +
              m_->AuGuPenalty(st2b, plb) + dp_.dp[piv + 1][en - 2][DP_U] + outer_coax;
          Expand(energy, {st + 2, piv, DP_P}, {piv + 1, en - 2, DP_U},
              {st + 2, CTD_LCOAX_WITH_PREV}, {en, CTD_LCOAX_WITH_NEXT});

          // (.   (   ).) Right outer coax
          energy = base_and_branch + dp_.dp[st + 2][piv][DP_U] + m_->multiloop_b +
              m_->AuGuPenalty(prb, en2b) + dp_.dp[piv + 1][en - 2][DP_P] + outer_coax;
          Expand(energy, {st + 2, piv, DP_U}, {piv + 1, en - 2, DP_P}, {piv + 1, CTD_RC_WITH_NEXT},
              {en, CTD_RC_WITH_PREV});

          // (.(   ).   ) Left inner coax
          energy = base_and_branch + dp_.dp[st + 2][piv - 1][DP_P] + m_->multiloop_b +
              m_->AuGuPenalty(st2b, pl1b) + dp_.dp[piv + 1][en - 1][DP_U] +
              m_->MismatchCoaxial(pl1b, plb, st1b, st2b) + m_->pf.Unpaired(st + 1) +
              m_->pf.Unpaired(piv);
          Expand(energy, {st + 2, piv - 1, DP_P}, {piv + 1, en - 1, DP_U},
              {st + 2, CTD_RC_WITH_PREV}, {en, CTD_RC_WITH_NEXT});

          // (   .(   ).) Right inner coax
          energy = base_and_branch + dp_.dp[st + 1][piv][DP_U] + m_->multiloop_b +
              m_->AuGuPenalty(pr1b, en2b) + dp_.dp[piv + 2][en - 2][DP_P] +
              m_->MismatchCoaxial(en2b, en1b, prb, pr1b) + m_->pf.Unpaired(piv + 1) +
              m_->pf.Unpaired(en - 1);
          Expand(energy, {st + 1, piv, DP_U}, {piv + 2, en - 2, DP_P},
              {piv + 2, CTD_LCOAX_WITH_NEXT}, {en, CTD_LCOAX_WITH_PREV});

          // ((   )   ) Left flush coax
          energy = base_and_branch + dp_.dp[st + 1][piv][DP_P] + m_->multiloop_b +
              m_->AuGuPenalty(st1b, plb) + dp_.dp[piv + 1][en - 1][DP_U] +
              m_->stack[stb][st1b][plb][enb];
          Expand(energy, {st + 1, piv, DP_P}, {piv + 1, en - 1, DP_U},
              {st + 1, CTD_FCOAX_WITH_PREV}, {en, CTD_FCOAX_WITH_NEXT});

          // (   (   )) Right flush coax
          energy = base_and_branch + dp_.dp[st + 1][piv][DP_U] + m_->multiloop_b +
              m_->AuGuPenalty(prb, en1b) + dp_.dp[piv + 1][en - 1][DP_P] +
              m_->stack[stb][prb][en1b][enb];
          Expand(energy, {st + 1, piv, DP_U}, {piv + 1, en - 1, DP_P},
              {piv + 1, CTD_FCOAX_WITH_NEXT}, {en, CTD_FCOAX_WITH_PREV});
        }
      }
    } else {
      // Left unpaired. Either DP_U or DP_U2.
      if (st + 1 < en && (a == DP_U || a == DP_U2)) {
        energy = base_energy + m_->pf.Unpaired(st) + dp_.dp[st + 1][en][a];
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
        auto base00 = dp_.dp[st][piv][DP_P] + m_->AuGuPenalty(stb, pb) + m_->multiloop_b;
        auto base01 = dp_.dp[st][piv - 1][DP_P] + m_->AuGuPenalty(stb, pl1b) + m_->multiloop_b;
        auto base10 = dp_.dp[st + 1][piv][DP_P] + m_->AuGuPenalty(st1b, pb) + m_->multiloop_b;
        auto base11 = dp_.dp[st + 1][piv - 1][DP_P] + m_->AuGuPenalty(st1b, pl1b) + m_->multiloop_b;

        const auto right_paired = dp_.dp[piv + 1][en][DP_U];
        // This is only usable if a != DP_U2 since this leaves everything unpaired.
        const auto right_unpaired = m_->pf.UnpairedCum(piv + 1, en);

        // Check a == U_RC:
        // (   )<.( ** ). > Right coax backward
        if (m_->cfg().UseCoaxialStacking() && a == DP_U_RC) {
          energy = base_energy + base11 + m_->MismatchCoaxial(pl1b, pb, stb, st1b) +
              m_->pf.Unpaired(st) + m_->pf.Unpaired(piv);
          // Our ctds will have already been set by now.
          Expand(energy + right_unpaired, {st + 1, piv - 1, DP_P});
          Expand(energy + right_paired, {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U});
        }

        // DP_U_RC is only the above case.
        if (a == DP_U_RC) continue;

        // From here on, a must be U, U2, U_WC, or U_GU.

        // (   )<   > - U, U2, U_WC?, U_GU?
        energy = base_energy + base00;
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
          Expand(energy + right_unpaired, {st, piv, DP_P}, {st, val_ctd});
          Expand(energy + right_paired, {st, piv, DP_P}, {piv + 1, en, DP_U}, {st, val_ctd});
        }

        if (a == DP_U2)
          Expand(energy + right_paired, {st, piv, DP_P}, {piv + 1, en, DP_U}, {st, val_ctd});

        // Make sure we don't form any branches that are not the right type of pair.
        if ((a == DP_U_WC && IsWcPair(stb, pb)) || (a == DP_U_GU && IsGuPair(stb, pb))) {
          Expand(energy + right_unpaired, {st, piv, DP_P});
          Expand(energy + right_paired, {st, piv, DP_P}, {piv + 1, en, DP_U});
        }

        // The rest of the cases are for U and U2.
        if (a != DP_U && a != DP_U2) continue;

        if (m_->cfg().UseDangleMismatch()) {
          // (   )3<   > 3' - U, U2
          energy = base_energy + base01 + m_->dangle3[pl1b][pb][stb] + m_->pf.Unpaired(piv);
          // Can only let the rest be unpaired if we only need one branch, i.e. DP_U not DP_U2.
          if (a == DP_U) Expand(energy + right_unpaired, {st, piv - 1, DP_P}, {st, CTD_3_DANGLE});
          Expand(
              energy + right_paired, {st, piv - 1, DP_P}, {piv + 1, en, DP_U}, {st, CTD_3_DANGLE});

          // 5(   )<   > 5' - U, U2
          energy = base_energy + base10 + m_->dangle5[pb][stb][st1b] + m_->pf.Unpaired(st);
          if (a == DP_U)
            Expand(energy + right_unpaired, {st + 1, piv, DP_P}, {st + 1, CTD_5_DANGLE});
          Expand(energy + right_paired, {st + 1, piv, DP_P}, {piv + 1, en, DP_U},
              {st + 1, CTD_5_DANGLE});

          // .(   ).<   > Terminal mismatch - U, U2
          energy = base_energy + base11 + m_->terminal[pl1b][pb][stb][st1b] + m_->pf.Unpaired(st) +
              m_->pf.Unpaired(piv);
          if (a == DP_U)
            Expand(energy + right_unpaired, {st + 1, piv - 1, DP_P}, {st + 1, CTD_MISMATCH});
          Expand(energy + right_paired, {st + 1, piv - 1, DP_P}, {piv + 1, en, DP_U},
              {st + 1, CTD_MISMATCH});
        }

        if (m_->cfg().UseCoaxialStacking()) {
          // .(   ).<(   ) > Left coax - U, U2
          energy = base_energy + base11 + m_->MismatchCoaxial(pl1b, pb, stb, st1b) +
              m_->pf.Unpaired(st) + m_->pf.Unpaired(piv);
          Expand(energy + dp_.dp[piv + 1][en][DP_U_WC], {st + 1, piv - 1, DP_P},
              {piv + 1, en, DP_U_WC}, {st + 1, CTD_LCOAX_WITH_NEXT},
              {piv + 1, CTD_LCOAX_WITH_PREV});
          Expand(energy + dp_.dp[piv + 1][en][DP_U_GU], {st + 1, piv - 1, DP_P},
              {piv + 1, en, DP_U_GU}, {st + 1, CTD_LCOAX_WITH_NEXT},
              {piv + 1, CTD_LCOAX_WITH_PREV});

          // (   )(<   ) > Flush coax - U, U2
          energy = base_energy + base01 + m_->stack[pl1b][pb][WcPair(pb)][stb] +
              dp_.dp[piv][en][DP_U_WC];
          Expand(energy, {st, piv - 1, DP_P}, {piv, en, DP_U_WC}, {st, CTD_FCOAX_WITH_NEXT},
              {piv, CTD_FCOAX_WITH_PREV});

          if (IsGu(pb)) {
            energy = base_energy + base01 + m_->stack[pl1b][pb][GuPair(pb)][stb] +
                dp_.dp[piv][en][DP_U_GU];
            Expand(energy, {st, piv - 1, DP_P}, {piv, en, DP_U_GU}, {st, CTD_FCOAX_WITH_NEXT},
                {piv, CTD_FCOAX_WITH_PREV});
          }

          if (piv < en - 1) {
            // (   )<.(   ). > Right coax forward - U, U2
            energy = base_energy + base00 + dp_.dp[piv + 1][en][DP_U_RC];
            Expand(energy, {st, piv, DP_P}, {piv + 1, en, DP_U_RC}, {st, CTD_RC_WITH_NEXT},
                {piv + 2, CTD_RC_WITH_PREV});
          }
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

}  // namespace mrna::md::base
