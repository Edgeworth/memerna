// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "backends/baseopt/energy/model.h"
#include "backends/baseopt/energy/precomp.h"
#include "backends/common/base/dp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/array.h"
#include "util/error.h"

namespace mrna::md::base::opt {

void MfeLyngsoSparseOpt(const Primary& r, const Model::Ptr& m, DpState& state) {
  static_assert(
      HAIRPIN_MIN_SZ >= 3, "Minimum hairpin size >= 3 is relied upon in some expressions.");

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false, true},
      .ctd{erg::EnergyCfg::Ctd::ALL},
  };
  support.VerifySupported(funcname(), m->cfg());

  spdlog::debug("baseopt {} with cfg {}", funcname(), m->cfg());

  const int N = static_cast<int>(r.size());
  const Precomp pc(Primary(r), m);
  state.dp = DpArray(r.size() + 1, MAX_E);
  auto& dp = state.dp;

  // See ComputeTables2 for comments - it is mostly the same.
  std::vector<std::vector<Cand>> p_cand_en[CAND_EN_SIZE];
  for (auto& i : p_cand_en) i.resize(r.size());
  std::vector<Cand> cand_st[CAND_SIZE];
  Array2D1S<Energy, TWOLOOP_MAX_SZ + 1> lyngso(r.size(), MAX_E);
  for (int st = N - 1; st >= 0; --st) {
    for (auto& i : cand_st) i.clear();
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const Base stb = r[st];
      const Base st1b = r[st + 1];
      const Base st2b = r[st + 2];
      const Base enb = r[en];
      const Base en1b = r[en - 1];
      const Base en2b = r[en - 2];
      Energy mins[] = {MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(sizeof(mins) / sizeof(mins[0]) == DP_SIZE, "array wrong size");
      const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);

      for (int l = 0; l <= max_inter; ++l) {
        // Don't add asymmetry here
        if (l >= 2)
          lyngso[st][en][l] = std::min(lyngso[st][en][l],
              lyngso[st + 1][en - 1][l - 2] - m->internal_init[l - 2] + m->internal_init[l]);

        // Add asymmetry here, on left and right
        auto val = std::min(l * m->internal_asym, NINIO_MAX_ASYM) + m->internal_init[l];
        lyngso[st][en][l] = std::min(lyngso[st][en][l],
            m->InternalLoopAuGuPenalty(r[st + l + 1], en1b) + m->AuGuPenalty(r[st + l + 1], en1b) +
                m->internal_other_mismatch[en1b][enb][r[st + l]][r[st + l + 1]] + val +
                dp[st + l + 1][en - 1][DP_P]);
        lyngso[st][en][l] = std::min(lyngso[st][en][l],
            m->InternalLoopAuGuPenalty(st1b, r[en - l - 1]) + m->AuGuPenalty(st1b, r[en - l - 1]) +
                m->internal_other_mismatch[r[en - l - 1]][r[en - l]][stb][st1b] + val +
                dp[st + 1][en - l - 1][DP_P]);
      }

      if (m->CanPair(r, st, en)) {
        // Stacking
        mins[DP_P] =
            std::min(mins[DP_P], m->stack[stb][st1b][en1b][enb] + dp[st + 1][en - 1][DP_P]);
        // Bulge
        for (int isz = 1; isz <= max_inter; ++isz) {
          mins[DP_P] = std::min(mins[DP_P],
              m->Bulge(r, st, en, st + 1 + isz, en - 1) + dp[st + 1 + isz][en - 1][DP_P]);
          mins[DP_P] = std::min(mins[DP_P],
              m->Bulge(r, st, en, st + 1, en - 1 - isz) + dp[st + 1][en - 1 - isz][DP_P]);
        }

        // Ax1 internal loops. Make sure to skip 0x1, 1x1, 2x1, and 1x2 loops, since they have
        // special energies.
        static_assert(Model::INITIATION_CACHE_SZ > TWOLOOP_MAX_SZ,
            "need initiation cached up to TWOLOOP_MAX_SZ");
        auto base_internal_loop = m->InternalLoopAuGuPenalty(stb, enb) + m->AuGuPenalty(stb, enb);
        for (int isz = 4; isz <= max_inter; ++isz) {
          auto val = base_internal_loop + m->internal_init[isz] +
              std::min((isz - 2) * m->internal_asym, NINIO_MAX_ASYM);
          mins[DP_P] = std::min(mins[DP_P],
              val + m->InternalLoopAuGuPenalty(r[st + isz], en2b) +
                  m->AuGuPenalty(r[st + isz], en2b) + dp[st + isz][en - 2][DP_P]);
          mins[DP_P] = std::min(mins[DP_P],
              val + m->InternalLoopAuGuPenalty(st2b, r[en - isz]) +
                  m->AuGuPenalty(st2b, r[en - isz]) + dp[st + 2][en - isz][DP_P]);
        }

        // Internal loop cases. Since we require HAIRPIN_MIN_SZ >= 3 and initialise arr to MAX_E, we
        // don't need ifs here.
        mins[DP_P] = std::min(mins[DP_P],
            m->internal_1x1[stb][st1b][st2b][en2b][en1b][enb] + m->AuGuPenalty(stb, enb) +
                m->AuGuPenalty(st2b, en2b) + dp[st + 2][en - 2][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            m->internal_1x2[stb][st1b][st2b][r[en - 3]][en2b][en1b][enb] +
                m->AuGuPenalty(stb, enb) + m->AuGuPenalty(st2b, r[en - 3]) +
                dp[st + 2][en - 3][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            m->internal_1x2[en2b][en1b][enb][stb][st1b][st2b][r[st + 3]] +
                m->AuGuPenalty(en2b, r[st + 3]) + m->AuGuPenalty(enb, stb) +
                dp[st + 3][en - 2][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            m->internal_2x2[stb][st1b][st2b][r[st + 3]][r[en - 3]][en2b][en1b][enb] +
                m->AuGuPenalty(stb, enb) + m->AuGuPenalty(r[st + 3], r[en - 3]) +
                dp[st + 3][en - 3][DP_P]);

        // 2x3 and 3x2 loops
        const auto two_by_three = base_internal_loop + m->internal_init[5] +
            std::min(m->internal_asym, NINIO_MAX_ASYM) +
            m->internal_2x3_mismatch[stb][st1b][en1b][enb];
        mins[DP_P] = std::min(mins[DP_P],
            two_by_three + m->InternalLoopAuGuPenalty(r[st + 3], r[en - 4]) +
                m->AuGuPenalty(r[st + 3], r[en - 4]) +
                m->internal_2x3_mismatch[r[en - 4]][r[en - 3]][st2b][r[st + 3]] +
                dp[st + 3][en - 4][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            two_by_three + m->InternalLoopAuGuPenalty(r[st + 4], r[en - 3]) +
                m->AuGuPenalty(r[st + 4], r[en - 3]) +
                m->internal_2x3_mismatch[r[en - 3]][r[en - 2]][r[st + 3]][r[st + 4]] +
                dp[st + 4][en - 3][DP_P]);

        // For the rest of the loops we need to apply the "other" type mismatches.
        base_internal_loop += m->internal_other_mismatch[stb][st1b][en1b][enb];

        // Lyngso for the rest.
        for (int l = 6; l <= max_inter; ++l)
          mins[DP_P] = std::min(mins[DP_P],
              lyngso[st + 2][en - 2][l - 4] - m->internal_init[l - 4] + m->internal_init[l] +
                  base_internal_loop);

        // Hairpin loops.
        mins[DP_P] = std::min(mins[DP_P], pc.Hairpin(st, en));

        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost = pc.augubranch[stb][enb] + m->multiloop_hack_a;

        // (<   ><   >)
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + dp[st + 1][en - 1][DP_U2]);

        // (3<   ><   >) 3'
        mins[DP_P] = std::min(
            mins[DP_P], base_branch_cost + dp[st + 2][en - 1][DP_U2] + m->dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        mins[DP_P] = std::min(
            mins[DP_P], base_branch_cost + dp[st + 1][en - 2][DP_U2] + m->dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        mins[DP_P] = std::min(mins[DP_P],
            base_branch_cost + dp[st + 2][en - 2][DP_U2] + m->terminal[stb][st1b][en1b][enb]);

        // (.(   ).   ) Left inner coax
        for (auto cand : cand_st[CAND_P_LIC])
          mins[DP_P] =
              std::min(mins[DP_P], base_branch_cost + cand.energy + dp[cand.idx][en - 1][DP_U]);
        // (.(   )   .) Left outer coax
        const auto outer_coax = m->MismatchCoaxial(stb, st1b, en1b, enb);
        for (auto cand : cand_st[CAND_P_LOC])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_mismatch_coax + outer_coax +
                  dp[cand.idx][en - 2][DP_U]);
        // ((   )   ) Left flush coax
        for (auto cand : cand_st[CAND_P_LFC])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_flush_coax +
                  m->stack[stb][st1b][r[cand.idx]][enb] + dp[cand.idx + 1][en - 1][DP_U]);
        // (   .(   ).) Right inner coax
        for (auto cand : p_cand_en[CAND_EN_P_RIC][en])
          mins[DP_P] =
              std::min(mins[DP_P], base_branch_cost + cand.energy + dp[st + 1][cand.idx][DP_U]);
        // (.   (   ).) Right outer coax
        for (auto cand : p_cand_en[CAND_EN_P_ROC][en])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_mismatch_coax + outer_coax +
                  dp[st + 2][cand.idx][DP_U]);
        // (   (   )) Right flush coax
        for (auto cand : p_cand_en[CAND_EN_P_RFC][en])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_flush_coax +
                  m->stack[stb][r[cand.idx]][en1b][enb] + dp[st + 1][cand.idx - 1][DP_U]);

        dp[st][en][DP_P] = mins[DP_P];
      }
      // Update unpaired.
      // Choose `st` to be unpaired.
      if (st + 1 < en) {
        mins[DP_U] = std::min(mins[DP_U], dp[st + 1][en][DP_U]);
        mins[DP_U2] = std::min(mins[DP_U2], dp[st + 1][en][DP_U2]);
      }
      for (auto cand : cand_st[CAND_U]) {
        mins[DP_U] = std::min(mins[DP_U], cand.energy + std::min(dp[cand.idx][en][DP_U], ZERO_E));
        mins[DP_U2] = std::min(mins[DP_U2], cand.energy + dp[cand.idx][en][DP_U]);
      }

      for (auto cand : cand_st[CAND_U_LC]) {
        // |.(   ).|<(   ) > Left coaxial stack
        const auto val =
            cand.energy + std::min(dp[cand.idx][en][DP_U_WC], dp[cand.idx][en][DP_U_GU]);
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_RC_FWD]) {
        // |(   )|<.(   ). > Right coax forward
        const auto val = cand.energy - pc.min_mismatch_coax + dp[cand.idx][en][DP_U_RC];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_LFC_WC]) {
        // (   )(<   ) > Flush coax - U
        const auto val = cand.energy + dp[cand.idx][en][DP_U_WC];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_LFC_GU]) {
        // (   )(<   ) > Flush coax - U
        const auto val = cand.energy + dp[cand.idx][en][DP_U_GU];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }

      for (auto cand : cand_st[CAND_U_WC])
        mins[DP_U_WC] =
            std::min(mins[DP_U_WC], cand.energy + std::min(dp[cand.idx][en][DP_U], ZERO_E));
      for (auto cand : cand_st[CAND_U_GU])
        mins[DP_U_GU] =
            std::min(mins[DP_U_GU], cand.energy + std::min(dp[cand.idx][en][DP_U], ZERO_E));
      for (auto cand : cand_st[CAND_U_RC]) {
        // (   )<.( * ). > Right coax backward
        mins[DP_U_RC] =
            std::min(mins[DP_U_RC], cand.energy + std::min(dp[cand.idx][en][DP_U], ZERO_E));
      }

      // Set these so we can use sparse folding.
      dp[st][en][DP_U] = mins[DP_U];
      dp[st][en][DP_U2] = mins[DP_U2];
      dp[st][en][DP_U_WC] = mins[DP_U_WC];
      dp[st][en][DP_U_GU] = mins[DP_U_GU];
      dp[st][en][DP_U_RC] = mins[DP_U_RC];

      // Now build the candidates arrays based off the current area [st, en].
      // In general, the idea is to see if there exists something we could replace a structure
      // with that is as good. e.g. we could replace a (   )3' with the equivalent U since we know
      // the energy for (   )3' that is, it is self contained. If replacing it with U[st][en] is
      // better, then we do not need to consider (...)3' when computing a larger U. In some cases
      // we use the minimum possible energy if we don't know the energy exactly for a structure
      // (e.g. RC). These orderings are useful to remember: U <= U_WC, U_GU, U2
      Energy cand_st_u = MAX_E;

      // Unpaired cases. These store the best pairs u_cand that begin at st.
      // begin means that the whole interaction starts at st. e.g. .(   ). starts one before the
      // paren.
      // (   ) - Normal - U, U2
      const auto normal_base = dp[st][en][DP_P] + pc.augubranch[stb][enb];
      if (normal_base < dp[st][en][DP_U] && normal_base < cand_st_u) cand_st_u = normal_base;

      // For U_GU and U_WC, they can't be replaced with DP_U, so we need to compare them to
      // something they can be
      // replaced with, i.e. themselves.
      if (IsGuPair(stb, enb)) {
        if (normal_base < CAP_E && normal_base < dp[st][en][DP_U_GU] &&
            (cand_st[CAND_U_GU].empty() || normal_base < cand_st[CAND_U_GU].back().energy))
          cand_st[CAND_U_GU].push_back({normal_base, en + 1});
        // Base case.
        dp[st][en][DP_U_GU] = std::min(dp[st][en][DP_U_GU], normal_base);
      } else {
        if (normal_base < CAP_E && normal_base < dp[st][en][DP_U_WC] &&
            (cand_st[CAND_U_WC].empty() || normal_base < cand_st[CAND_U_WC].back().energy))
          cand_st[CAND_U_WC].push_back({normal_base, en + 1});
        // Base case.
        dp[st][en][DP_U_WC] = std::min(dp[st][en][DP_U_WC], normal_base);
      }

      // Can only merge candidate lists for monotonicity if
      // the right part of the pivot is the same (from the same array).
      // Can only apply monotonicity optimisation to ones ending with min(U, 0).
      // (   ). - 3' - U, U2
      const auto dangle3_base =
          dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + m->dangle3[en1b][enb][stb];
      if (dangle3_base < dp[st][en][DP_U] && dangle3_base < cand_st_u) cand_st_u = dangle3_base;
      // .(   ) - 5' - U, U2
      const auto dangle5_base =
          dp[st + 1][en][DP_P] + pc.augubranch[st1b][enb] + m->dangle5[enb][stb][st1b];
      if (dangle5_base < dp[st][en][DP_U] && dangle5_base < cand_st_u) cand_st_u = dangle5_base;
      // .(   ). - Terminal mismatch - U, U2
      const auto terminal_base =
          dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] + m->terminal[en1b][enb][stb][st1b];
      if (terminal_base < dp[st][en][DP_U] && terminal_base < cand_st_u) cand_st_u = terminal_base;

      // Add potentials to the candidate lists.
      if (cand_st_u < CAP_E &&
          (cand_st[CAND_U].empty() || cand_st_u < cand_st[CAND_U].back().energy))
        cand_st[CAND_U].push_back({cand_st_u, en + 1});

      // .(   ).<(   ) > - Left coax - U, U2
      const auto lcoax_base = dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] +
          m->MismatchCoaxial(en1b, enb, stb, st1b);
      if (lcoax_base < CAP_E && lcoax_base < dp[st][en][DP_U])
        cand_st[CAND_U_LC].push_back({lcoax_base, en + 1});

      // (   )<.(   ). > Right coax forward - U, U2
      const auto rcoaxf_base = dp[st][en][DP_P] + pc.augubranch[stb][enb] + pc.min_mismatch_coax;
      if (rcoaxf_base < CAP_E && rcoaxf_base < dp[st][en][DP_U])
        cand_st[CAND_U_RC_FWD].push_back({rcoaxf_base, en + 1});

      // (   )<.( * ). > Right coax backward - RC
      // Again, we can't replace RC with U, we'd have to replace it with RC, so compare to
      // itself.
      const auto rcoaxb_base = dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] +
          m->MismatchCoaxial(en1b, enb, stb, st1b);
      if (rcoaxb_base < CAP_E && rcoaxb_base < dp[st][en][DP_U_RC] &&
          (cand_st[CAND_U_RC].empty() || rcoaxb_base < cand_st[CAND_U_RC].back().energy))
        cand_st[CAND_U_RC].push_back({rcoaxb_base, en + 1});
      // Base case.
      dp[st][en][DP_U_RC] = std::min(dp[st][en][DP_U_RC], rcoaxb_base);

      // (   )(<   ) > Flush coax - U, U2
      const auto wc_flush_base =
          dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + m->stack[en1b][enb][WcPair(enb)][stb];
      if (wc_flush_base < CAP_E && wc_flush_base < dp[st][en - 1][DP_U])
        cand_st[CAND_U_LFC_WC].push_back({wc_flush_base, en});
      if (IsGu(enb)) {
        const auto gu_flush_base =
            dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + m->stack[en1b][enb][GuPair(enb)][stb];
        if (gu_flush_base < CAP_E && gu_flush_base < dp[st][en - 1][DP_U])
          cand_st[CAND_U_LFC_GU].push_back({gu_flush_base, en});
      }

      // Paired cases
      // (.(   )   .) Left outer coax - P
      // Since we assumed the minimum energy coax stack and made this structure self contained,
      // we could potentially replace it with U[st + 1][en].
      const auto plocoax_base =
          dp[st + 2][en][DP_P] + pc.augubranch[st2b][enb] + pc.min_mismatch_coax;
      if (plocoax_base < CAP_E && plocoax_base < dp[st + 1][en][DP_U])
        cand_st[CAND_P_LOC].push_back({plocoax_base, en + 1});
      // (.   (   ).) Right outer coax
      const auto procoax_base =
          dp[st][en - 2][DP_P] + pc.augubranch[stb][en2b] + pc.min_mismatch_coax;
      if (procoax_base < CAP_E && procoax_base < dp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_ROC][en].push_back({procoax_base, st - 1});
      // (.(   ).   ) Left inner coax
      const auto plrcoax_base = dp[st + 2][en - 1][DP_P] + pc.augubranch[st2b][en1b] +
          m->MismatchCoaxial(en1b, enb, st1b, st2b);
      if (plrcoax_base < CAP_E && plrcoax_base < dp[st + 1][en][DP_U])
        cand_st[CAND_P_LIC].push_back({plrcoax_base, en + 1});
      // (   .(   ).) Right inner coax
      const auto prlcoax_base = dp[st + 1][en - 2][DP_P] + pc.augubranch[st1b][en2b] +
          m->MismatchCoaxial(en2b, en1b, stb, st1b);
      if (prlcoax_base < CAP_E && prlcoax_base < dp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_RIC][en].push_back({prlcoax_base, st - 1});
      // ((   )   ) Left flush coax
      const auto plfcoax_base = dp[st + 1][en][DP_P] + pc.augubranch[st1b][enb] + pc.min_flush_coax;
      if (plfcoax_base < CAP_E && plfcoax_base < dp[st + 1][en][DP_U])
        cand_st[CAND_P_LFC].push_back({plfcoax_base, en});
      // (   (   )) Right flush coax
      const auto prfcoax_base = dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + pc.min_flush_coax;
      if (prfcoax_base < CAP_E && prfcoax_base < dp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_RFC][en].push_back({prfcoax_base, st});

      // Base cases.
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], cand_st_u);
      // Note we don't include the stacking here since they can't be base cases for U.
    }
  }
}

}  // namespace mrna::md::base::opt
