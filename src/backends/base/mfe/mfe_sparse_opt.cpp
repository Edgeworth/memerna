// Copyright 2016 E.
#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/base/energy/precomp.h"
#include "backends/common/base/dp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/error.h"

namespace mrna::md::base {

void MfeSparseOpt(const Primary& r, const Model::Ptr& m, DpState& state) {
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  static thread_local const erg::EnergyCfgSupport support{
      .lonely_pairs{erg::EnergyCfg::LonelyPairs::HEURISTIC, erg::EnergyCfg::LonelyPairs::ON},
      .bulge_states{false, true},
      .ctd{erg::EnergyCfg::Ctd::ALL, erg::EnergyCfg::Ctd::NO_COAX, erg::EnergyCfg::Ctd::D2,
          erg::EnergyCfg::Ctd::NONE},
  };
  support.VerifySupported(funcname(), m->cfg());
  m->pf.Verify(r);

  spdlog::debug("base {} with cfg {}", funcname(), m->cfg());

  const int N = static_cast<int>(r.size());
  const Precomp pc(Primary(r), m);
  state.dp = DpArray(r.size() + 1, MAX_E);
  auto& dp = state.dp;
  std::vector<std::vector<Cand>> p_cand_en[CAND_EN_SIZE];
  for (auto& i : p_cand_en) i.resize(r.size());
  std::vector<Cand> cand_st[CAND_SIZE];
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

      if (m->CanPair(r, st, en)) {
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        mins[DP_P] = std::min(mins[DP_P],
            m->stack[stb][st1b][en1b][enb] + m->pf.Paired(st, en) + dp[st + 1][en - 1][DP_P]);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (dp[ist][ien][DP_P] < mins[DP_P] - pc.min_twoloop_not_stack - pc.sum_neg_pf)
              mins[DP_P] = std::min(mins[DP_P], pc.TwoLoop(st, en, ist, ien) + dp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        mins[DP_P] = std::min(mins[DP_P], pc.Hairpin(st, en));

        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost =
            pc.augubranch[stb][enb] + m->pf.Paired(st, en) + m->multiloop_a;

        // (<   ><   >)
        auto val = base_branch_cost + dp[st + 1][en - 1][DP_U2];

        if (m->cfg().UseD2()) {
          // D2 can overlap terminal mismatches with anything.
          // (<   ><   >) Terminal mismatch
          val += m->terminal[stb][st1b][en1b][enb];
        }

        mins[DP_P] = std::min(mins[DP_P], val);

        if (m->cfg().UseDangleMismatch()) {
          // (3<   ><   >) 3'
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + dp[st + 2][en - 1][DP_U2] + m->dangle3[stb][st1b][enb] +
                  m->pf.Unpaired(st + 1));
          // (<   ><   >5) 5'
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + dp[st + 1][en - 2][DP_U2] + m->dangle5[stb][en1b][enb] +
                  m->pf.Unpaired(en - 1));
          // (.<   ><   >.) Terminal mismatch
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + dp[st + 2][en - 2][DP_U2] + m->terminal[stb][st1b][en1b][enb] +
                  m->pf.Unpaired(st + 1) + m->pf.Unpaired(en - 1));
        }

        if (m->cfg().UseCoaxialStacking()) {
          // (.(   ).   ) Left inner coax
          for (auto cand : cand_st[CAND_P_LIC])
            mins[DP_P] =
                std::min(mins[DP_P], base_branch_cost + cand.energy + dp[cand.idx][en - 1][DP_U]);
          // (.(   )   .) Left outer coax
          const auto outer_coax = m->MismatchCoaxial(stb, st1b, en1b, enb);
          for (auto cand : cand_st[CAND_P_LOC])
            mins[DP_P] = std::min(mins[DP_P],
                base_branch_cost + cand.energy - pc.min_mismatch_coax - pc.min_pf_unpaired +
                    outer_coax + m->pf.Unpaired(en - 1) + dp[cand.idx][en - 2][DP_U]);
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
                base_branch_cost + cand.energy - pc.min_mismatch_coax - pc.min_pf_unpaired +
                    outer_coax + m->pf.Unpaired(st + 1) + dp[st + 2][cand.idx][DP_U]);
          // (   (   )) Right flush coax
          for (auto cand : p_cand_en[CAND_EN_P_RFC][en])
            mins[DP_P] = std::min(mins[DP_P],
                base_branch_cost + cand.energy - pc.min_flush_coax +
                    m->stack[stb][r[cand.idx]][en1b][enb] + dp[st + 1][cand.idx - 1][DP_U]);
        }

        dp[st][en][DP_P] = mins[DP_P];
      }
      // Update unpaired.
      // Choose `st` to be unpaired.
      if (st + 1 < en) {
        mins[DP_U] = std::min(mins[DP_U], dp[st + 1][en][DP_U] + m->pf.Unpaired(st));
        mins[DP_U2] = std::min(mins[DP_U2], dp[st + 1][en][DP_U2] + m->pf.Unpaired(st));
      }
      for (auto cand : cand_st[CAND_U]) {
        mins[DP_U] = std::min(mins[DP_U],
            cand.energy + std::min(dp[cand.idx][en][DP_U], m->pf.UnpairedCum(cand.idx, en)));
        mins[DP_U2] = std::min(mins[DP_U2], cand.energy + dp[cand.idx][en][DP_U]);
      }

      if (m->cfg().UseCoaxialStacking()) {
        for (auto cand : cand_st[CAND_U_LC]) {
          // |.(   ).|<(   ) > Left coaxial stack
          const auto val =
              cand.energy + std::min(dp[cand.idx][en][DP_U_WC], dp[cand.idx][en][DP_U_GU]);
          mins[DP_U] = std::min(mins[DP_U], val);
          mins[DP_U2] = std::min(mins[DP_U2], val);
        }
        for (auto cand : cand_st[CAND_U_RC_FWD]) {
          // |(   )|<.(   ). > Right coax forward
          const auto val = cand.energy - pc.min_mismatch_coax - pc.min_pf_unpaired * 2 +
              dp[cand.idx][en][DP_U_RC];
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
      }

      for (auto cand : cand_st[CAND_U_WC])
        mins[DP_U_WC] = std::min(mins[DP_U_WC],
            cand.energy + std::min(dp[cand.idx][en][DP_U], m->pf.UnpairedCum(cand.idx, en)));
      for (auto cand : cand_st[CAND_U_GU])
        mins[DP_U_GU] = std::min(mins[DP_U_GU],
            cand.energy + std::min(dp[cand.idx][en][DP_U], m->pf.UnpairedCum(cand.idx, en)));
      for (auto cand : cand_st[CAND_U_RC]) {
        // (   )<.( * ). > Right coax backward
        mins[DP_U_RC] = std::min(mins[DP_U_RC],
            cand.energy + std::min(dp[cand.idx][en][DP_U], m->pf.UnpairedCum(cand.idx, en)));
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
      auto normal_base = dp[st][en][DP_P] + pc.augubranch[stb][enb];

      if (m->cfg().UseD2()) {
        // Note that D2 can overlap with anything.
        if (st != 0 && en != N - 1) {
          // (   )<   > Terminal mismatch - U
          normal_base += m->terminal[enb][r[en + 1]][r[st - 1]][stb];
        } else if (en != N - 1) {
          // (   )<3   > 3' - U
          normal_base += m->dangle3[enb][r[en + 1]][stb];
        } else if (st != 0) {
          // 5(   )<   > 5' - U
          normal_base += m->dangle5[enb][r[st - 1]][stb];
        }
      }

      if (normal_base < dp[st][en][DP_U] && normal_base < cand_st_u) cand_st_u = normal_base;

      // For U_GU and U_WC, they can't be replaced with DP_U, so we need to compare them to
      // something they can be
      // replaced with, i.e. themselves.
      if (IsGuPair(stb, enb)) {
        const auto* cand_st_u_gu_last =
            cand_st[CAND_U_GU].empty() ? nullptr : &cand_st[CAND_U_GU].back();
        const Energy cand_st_u_gu_unpaired_cum =
            cand_st_u_gu_last ? m->pf.UnpairedCum(cand_st_u_gu_last->idx, en) : ZERO_E;
        if (normal_base < CAP_E && normal_base < dp[st][en][DP_U_GU] &&
            (!cand_st_u_gu_last ||
                normal_base < cand_st_u_gu_last->energy + cand_st_u_gu_unpaired_cum))
          cand_st[CAND_U_GU].push_back({normal_base, en + 1});
        // Base case.
        dp[st][en][DP_U_GU] = std::min(dp[st][en][DP_U_GU], normal_base);
      } else {
        const auto* cand_st_u_wc_last =
            cand_st[CAND_U_WC].empty() ? nullptr : &cand_st[CAND_U_WC].back();
        const Energy cand_st_u_wc_unpaired_cum =
            cand_st_u_wc_last ? m->pf.UnpairedCum(cand_st_u_wc_last->idx, en) : ZERO_E;
        if (normal_base < CAP_E && normal_base < dp[st][en][DP_U_WC] &&
            (!cand_st_u_wc_last ||
                normal_base < cand_st_u_wc_last->energy + cand_st_u_wc_unpaired_cum))
          cand_st[CAND_U_WC].push_back({normal_base, en + 1});
        // Base case.
        dp[st][en][DP_U_WC] = std::min(dp[st][en][DP_U_WC], normal_base);
      }

      // Can only merge candidate lists for monotonicity if
      // the right part of the pivot is the same (from the same array).
      // Can only apply monotonicity optimisation to ones ending with min(U, 0).
      if (m->cfg().UseDangleMismatch()) {
        // (   ). - 3' - U, U2
        const auto dangle3_base = dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] +
            m->dangle3[en1b][enb][stb] + m->pf.Unpaired(en);
        if (dangle3_base < dp[st][en][DP_U] && dangle3_base < cand_st_u) cand_st_u = dangle3_base;
        // .(   ) - 5' - U, U2
        const auto dangle5_base = dp[st + 1][en][DP_P] + pc.augubranch[st1b][enb] +
            m->dangle5[enb][stb][st1b] + m->pf.Unpaired(st);
        if (dangle5_base < dp[st][en][DP_U] && dangle5_base < cand_st_u) cand_st_u = dangle5_base;
        // .(   ). - Terminal mismatch - U, U2
        const auto terminal_base = dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] +
            m->terminal[en1b][enb][stb][st1b] + m->pf.Unpaired(st) + m->pf.Unpaired(en);
        if (terminal_base < dp[st][en][DP_U] && terminal_base < cand_st_u)
          cand_st_u = terminal_base;
      }

      // Add potentials to the candidate lists. Need to adjust by the unpaired pseudofree energy if
      // it exists to maintain the monotonicity property.
      const auto* cand_st_u_last = cand_st[CAND_U].empty() ? nullptr : &cand_st[CAND_U].back();
      const Energy cand_st_u_unpaired_cum =
          cand_st_u_last ? m->pf.UnpairedCum(cand_st_u_last->idx, en) : ZERO_E;
      if (cand_st_u < CAP_E &&
          (!cand_st_u_last || cand_st_u < cand_st_u_last->energy + cand_st_u_unpaired_cum))
        cand_st[CAND_U].push_back({cand_st_u, en + 1});

      if (m->cfg().UseCoaxialStacking()) {
        // .(   ).<(   ) > - Left coax - U, U2
        const auto lcoax_base = dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] +
            m->MismatchCoaxial(en1b, enb, stb, st1b) + m->pf.Unpaired(st) + m->pf.Unpaired(en);
        if (lcoax_base < CAP_E && lcoax_base < dp[st][en][DP_U])
          cand_st[CAND_U_LC].push_back({lcoax_base, en + 1});

        // (   )<.(   ). > Right coax forward - U, U2
        const auto rcoaxf_base = dp[st][en][DP_P] + pc.augubranch[stb][enb] + pc.min_mismatch_coax +
            pc.min_pf_unpaired * 2;
        if (rcoaxf_base < CAP_E && rcoaxf_base < dp[st][en][DP_U])
          cand_st[CAND_U_RC_FWD].push_back({rcoaxf_base, en + 1});

        // (   )<.( * ). > Right coax backward - RC
        // Again, we can't replace RC with U, we'd have to replace it with RC, so compare to
        // itself.
        const auto rcoaxb_base = dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] +
            m->MismatchCoaxial(en1b, enb, stb, st1b) + m->pf.Unpaired(st) + m->pf.Unpaired(en);
        const auto* cand_st_u_rc_last =
            cand_st[CAND_U_RC].empty() ? nullptr : &cand_st[CAND_U_RC].back();
        const Energy cand_st_u_rc_unpaired_cum =
            cand_st_u_rc_last ? m->pf.UnpairedCum(cand_st_u_rc_last->idx, en) : ZERO_E;
        if (rcoaxb_base < CAP_E && rcoaxb_base < dp[st][en][DP_U_RC] &&
            (!cand_st_u_rc_last ||
                rcoaxb_base < cand_st_u_rc_last->energy + cand_st_u_rc_unpaired_cum))
          cand_st[CAND_U_RC].push_back({rcoaxb_base, en + 1});
        // Base case.
        dp[st][en][DP_U_RC] = std::min(dp[st][en][DP_U_RC], rcoaxb_base);

        // (   )(<   ) > Flush coax - U, U2
        const auto wc_flush_base =
            dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + m->stack[en1b][enb][WcPair(enb)][stb];
        if (wc_flush_base < CAP_E && wc_flush_base < dp[st][en - 1][DP_U])
          cand_st[CAND_U_LFC_WC].push_back({wc_flush_base, en});
        if (IsGu(enb)) {
          const auto gu_flush_base = dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] +
              m->stack[en1b][enb][GuPair(enb)][stb];
          if (gu_flush_base < CAP_E && gu_flush_base < dp[st][en - 1][DP_U])
            cand_st[CAND_U_LFC_GU].push_back({gu_flush_base, en});
        }

        // Paired cases
        // (.(   )   .) Left outer coax - P
        // Since we assumed the minimum energy coax stack and made this structure self contained,
        // we could potentially replace it with U[st + 1][en].
        const auto plocoax_base = dp[st + 2][en][DP_P] + pc.augubranch[st2b][enb] +
            m->pf.Unpaired(st + 1) + pc.min_mismatch_coax + pc.min_pf_unpaired;
        if (plocoax_base < CAP_E && plocoax_base < dp[st + 1][en][DP_U])
          cand_st[CAND_P_LOC].push_back({plocoax_base, en + 1});
        // (.   (   ).) Right outer coax
        const auto procoax_base = dp[st][en - 2][DP_P] + pc.augubranch[stb][en2b] +
            m->pf.Unpaired(en - 1) + pc.min_mismatch_coax + pc.min_pf_unpaired;
        if (procoax_base < CAP_E && procoax_base < dp[st][en - 1][DP_U])
          p_cand_en[CAND_EN_P_ROC][en].push_back({procoax_base, st - 1});
        // (.(   ).   ) Left inner coax
        const auto plrcoax_base = dp[st + 2][en - 1][DP_P] + pc.augubranch[st2b][en1b] +
            m->MismatchCoaxial(en1b, enb, st1b, st2b) + m->pf.Unpaired(st + 1) + m->pf.Unpaired(en);
        if (plrcoax_base < CAP_E && plrcoax_base < dp[st + 1][en][DP_U])
          cand_st[CAND_P_LIC].push_back({plrcoax_base, en + 1});
        // (   .(   ).) Right inner coax
        const auto prlcoax_base = dp[st + 1][en - 2][DP_P] + pc.augubranch[st1b][en2b] +
            m->MismatchCoaxial(en2b, en1b, stb, st1b) + m->pf.Unpaired(st) + m->pf.Unpaired(en - 1);
        if (prlcoax_base < CAP_E && prlcoax_base < dp[st][en - 1][DP_U])
          p_cand_en[CAND_EN_P_RIC][en].push_back({prlcoax_base, st - 1});
        // ((   )   ) Left flush coax
        const auto plfcoax_base =
            dp[st + 1][en][DP_P] + pc.augubranch[st1b][enb] + pc.min_flush_coax;
        if (plfcoax_base < CAP_E && plfcoax_base < dp[st + 1][en][DP_U])
          cand_st[CAND_P_LFC].push_back({plfcoax_base, en});
        // (   (   )) Right flush coax
        const auto prfcoax_base =
            dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + pc.min_flush_coax;
        if (prfcoax_base < CAP_E && prfcoax_base < dp[st][en - 1][DP_U])
          p_cand_en[CAND_EN_P_RFC][en].push_back({prfcoax_base, st});
      }

      // Base cases.
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], cand_st_u);
      // Note we don't include the stacking here since they can't be base cases for U.
    }
  }
}

}  // namespace mrna::md::base
