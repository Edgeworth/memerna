// Copyright 2016 Eliot Courtney.
#include <algorithm>
#include <memory>
#include <vector>

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "compute/energy/precomp.h"
#include "compute/mfe/mfe.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/primary.h"
#include "util/array.h"

namespace mrna::mfe {

DpArray ComputeTables2(const Primary& r, const energy::EnergyModelPtr& em) {
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  const int N = static_cast<int>(r.size());
  const energy::Precomp pc(Primary(r), em);
  auto dp = DpArray(r.size() + 1, MAX_E);

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

      if (em->CanPair(r, st, en)) {
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        mins[DP_P] =
            std::min(mins[DP_P], em->stack[stb][st1b][en1b][enb] + dp[st + 1][en - 1][DP_P]);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (dp[ist][ien][DP_P] < mins[DP_P] - pc.min_twoloop_not_stack)
              mins[DP_P] = std::min(mins[DP_P], pc.TwoLoop(st, en, ist, ien) + dp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        mins[DP_P] = std::min(mins[DP_P], pc.Hairpin(st, en));

        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost = pc.augubranch[stb][enb] + em->multiloop_hack_a;

        // (<   ><   >)
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + dp[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        mins[DP_P] = std::min(
            mins[DP_P], base_branch_cost + dp[st + 2][en - 1][DP_U2] + em->dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        mins[DP_P] = std::min(
            mins[DP_P], base_branch_cost + dp[st + 1][en - 2][DP_U2] + em->dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        mins[DP_P] = std::min(mins[DP_P],
            base_branch_cost + dp[st + 2][en - 2][DP_U2] + em->terminal[stb][st1b][en1b][enb]);

        // (.(   ).   ) Left inner coax
        for (auto cand : cand_st[CAND_P_LIC])
          mins[DP_P] =
              std::min(mins[DP_P], base_branch_cost + cand.energy + dp[cand.idx][en - 1][DP_U]);
        // (.(   )   .) Left outer coax
        const auto outer_coax = em->MismatchCoaxial(stb, st1b, en1b, enb);
        for (auto cand : cand_st[CAND_P_LOC])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_mismatch_coax + outer_coax +
                  dp[cand.idx][en - 2][DP_U]);
        // ((   )   ) Left flush coax
        for (auto cand : cand_st[CAND_P_LFC])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_flush_coax +
                  em->stack[stb][st1b][r[cand.idx]][enb] + dp[cand.idx + 1][en - 1][DP_U]);
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
                  em->stack[stb][r[cand.idx]][en1b][enb] + dp[st + 1][cand.idx - 1][DP_U]);

        dp[st][en][DP_P] = mins[DP_P];
      }
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        mins[DP_U] = std::min(mins[DP_U], dp[st + 1][en][DP_U]);
        mins[DP_U2] = std::min(mins[DP_U2], dp[st + 1][en][DP_U2]);
      }
      for (auto cand : cand_st[CAND_U]) {
        mins[DP_U] = std::min(mins[DP_U], cand.energy + std::min(dp[cand.idx][en][DP_U], ZERO_E));
        mins[DP_U2] = std::min(mins[DP_U2], cand.energy + dp[cand.idx][en][DP_U]);
      }
      for (auto cand : cand_st[CAND_U_LC]) {
        const auto val =
            cand.energy + std::min(dp[cand.idx][en][DP_U_WC], dp[cand.idx][en][DP_U_GU]);
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_RC_FWD]) {
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
        auto val = cand.energy + dp[cand.idx][en][DP_U_GU];
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
      // In general, the idea is to see if there exists something we could replace a structure with
      // that is as good.
      // e.g. we could replace a (   )3' with the equivalent U since we know the energy for (   )3'
      // that is, it is self contained. If replacing it with U[st][en] is better, then we do not
      // need to consider (...)3' when computing a larger U. In some cases we use the minimum
      // possible energy if we don't know the energy exactly for a structure (e.g. RC). These
      // orderings are useful to remember: U <= U_WC, U_GU, U2
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
          dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + em->dangle3[en1b][enb][stb];
      if (dangle3_base < dp[st][en][DP_U] && dangle3_base < cand_st_u) cand_st_u = dangle3_base;
      // .(   ) - 5' - U, U2
      const auto dangle5_base =
          dp[st + 1][en][DP_P] + pc.augubranch[st1b][enb] + em->dangle5[enb][stb][st1b];
      if (dangle5_base < dp[st][en][DP_U] && dangle5_base < cand_st_u) cand_st_u = dangle5_base;
      // .(   ). - Terminal mismatch - U, U2
      const auto terminal_base =
          dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] + em->terminal[en1b][enb][stb][st1b];
      if (terminal_base < dp[st][en][DP_U] && terminal_base < cand_st_u) cand_st_u = terminal_base;

      // Add potentials to the candidate lists.
      if (cand_st_u < CAP_E &&
          (cand_st[CAND_U].empty() || cand_st_u < cand_st[CAND_U].back().energy))
        cand_st[CAND_U].push_back({cand_st_u, en + 1});

      // .(   ).<(   ) > - Left coax - U, U2
      const auto lcoax_base = dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] +
          em->MismatchCoaxial(en1b, enb, stb, st1b);
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
          em->MismatchCoaxial(en1b, enb, stb, st1b);
      if (rcoaxb_base < CAP_E && rcoaxb_base < dp[st][en][DP_U_RC] &&
          (cand_st[CAND_U_RC].empty() || rcoaxb_base < cand_st[CAND_U_RC].back().energy))
        cand_st[CAND_U_RC].push_back({rcoaxb_base, en + 1});
      // Base case.
      dp[st][en][DP_U_RC] = std::min(dp[st][en][DP_U_RC], rcoaxb_base);

      // (   )(<   ) > Flush coax - U, U2
      const auto wc_flush_base =
          dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + em->stack[en1b][enb][WcPair(enb)][stb];
      if (wc_flush_base < CAP_E && wc_flush_base < dp[st][en - 1][DP_U])
        cand_st[CAND_U_LFC_WC].push_back({wc_flush_base, en});
      if (IsGu(enb)) {
        const auto gu_flush_base = dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] +
            em->stack[en1b][enb][GuPair(enb)][stb];
        if (gu_flush_base < CAP_E && gu_flush_base < dp[st][en - 1][DP_U])
          cand_st[CAND_U_LFC_GU].push_back({gu_flush_base, en});
      }

      // Base cases.
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], normal_base);
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], dangle3_base);
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], dangle5_base);
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], terminal_base);
      // Note we don't include the stacking here since they can't be base cases for U.

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
          em->MismatchCoaxial(en1b, enb, st1b, st2b);
      if (plrcoax_base < CAP_E && plrcoax_base < dp[st + 1][en][DP_U])
        cand_st[CAND_P_LIC].push_back({plrcoax_base, en + 1});
      // (   .(   ).) Right inner coax
      const auto prlcoax_base = dp[st + 1][en - 2][DP_P] + pc.augubranch[st1b][en2b] +
          em->MismatchCoaxial(en2b, en1b, stb, st1b);
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
    }
  }

  return dp;
}

}  // namespace mrna::mfe
