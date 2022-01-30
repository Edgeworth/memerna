// Copyright 2016 E.
#include <algorithm>
#include <memory>
#include <vector>

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "compute/energy/precomp.h"
#include "compute/mfe/mfe.h"
#include "model/base.h"
#include "model/model.h"
#include "model/primary.h"
#include "util/array.h"

namespace mrna::mfe {

DpArray ComputeTables3(const Primary& r, const energy::EnergyModel& em) {
  static_assert(
      HAIRPIN_MIN_SZ >= 3, "Minimum hairpin size >= 3 is relied upon in some expressions.");

  const int N = static_cast<int>(r.size());
  const energy::Precomp pc(Primary(r), em);
  auto dp = DpArray(r.size() + 1, MAX_E);

  // See ComputeTables2 for comments - it is mostly the same.
  std::vector<std::vector<Cand>> p_cand_en[CAND_EN_SIZE];
  for (auto& i : p_cand_en) i.resize(r.size());
  std::vector<Cand> cand_st[CAND_SIZE];
  Array2D1S<Energy, TWOLOOP_MAX_SZ + 1> lyngso(r.size(), MAX_E);
  for (int st = N - 1; st >= 0; --st) {
    for (auto& i : cand_st) i.clear();
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const Base stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1],
                 en2b = r[en - 2];
      Energy mins[] = {MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(sizeof(mins) / sizeof(mins[0]) == DP_SIZE, "array wrong size");
      const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);

      for (int l = 0; l <= max_inter; ++l) {
        // Don't add asymmetry here
        if (l >= 2)
          lyngso[st][en][l] = std::min(lyngso[st][en][l],
              lyngso[st + 1][en - 1][l - 2] - em.internal_init[l - 2] + em.internal_init[l]);

        // Add asymmetry here, on left and right
        auto val = std::min(l * em.internal_asym, NINIO_MAX_ASYM) + em.internal_init[l];
        lyngso[st][en][l] = std::min(lyngso[st][en][l],
            em.InternalLoopAuGuPenalty(r[st + l + 1], en1b) +
                em.internal_other_mismatch[en1b][enb][r[st + l]][r[st + l + 1]] + val +
                dp[st + l + 1][en - 1][DP_P]);
        lyngso[st][en][l] = std::min(lyngso[st][en][l],
            em.InternalLoopAuGuPenalty(st1b, r[en - l - 1]) +
                em.internal_other_mismatch[r[en - l - 1]][r[en - l]][stb][st1b] + val +
                dp[st + 1][en - l - 1][DP_P]);
      }

      // TODO: check lonely pairs
      if (em.CanPair(r, st, en)) {
        // Stacking
        mins[DP_P] =
            std::min(mins[DP_P], em.stack[stb][st1b][en1b][enb] + dp[st + 1][en - 1][DP_P]);
        // Bulge
        for (int isz = 1; isz <= max_inter; ++isz) {
          mins[DP_P] = std::min(mins[DP_P],
              em.Bulge(r, st, en, st + 1 + isz, en - 1) + dp[st + 1 + isz][en - 1][DP_P]);
          mins[DP_P] = std::min(mins[DP_P],
              em.Bulge(r, st, en, st + 1, en - 1 - isz) + dp[st + 1][en - 1 - isz][DP_P]);
        }

        // Ax1 internal loops. Make sure to skip 0x1, 1x1, 2x1, and 1x2 loops, since they have
        // special energies.
        static_assert(energy::EnergyModel::INITIATION_CACHE_SZ > TWOLOOP_MAX_SZ,
            "need initiation cached up to TWOLOOP_MAX_SZ");
        auto base_internal_loop = em.InternalLoopAuGuPenalty(stb, enb);
        for (int isz = 4; isz <= max_inter; ++isz) {
          auto val = base_internal_loop + em.internal_init[isz] +
              std::min((isz - 2) * em.internal_asym, NINIO_MAX_ASYM);
          mins[DP_P] = std::min(mins[DP_P],
              val + em.InternalLoopAuGuPenalty(r[st + isz], en2b) + dp[st + isz][en - 2][DP_P]);
          mins[DP_P] = std::min(mins[DP_P],
              val + em.InternalLoopAuGuPenalty(st2b, r[en - isz]) + dp[st + 2][en - isz][DP_P]);
        }

        // Internal loop cases. Since we require HAIRPIN_MIN_SZ >= 3 and initialise arr to MAX_E, we
        // don't need ifs
        // here.
        mins[DP_P] = std::min(mins[DP_P],
            em.internal_1x1[stb][st1b][st2b][en2b][en1b][enb] + dp[st + 2][en - 2][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            em.internal_1x2[stb][st1b][st2b][r[en - 3]][en2b][en1b][enb] +
                dp[st + 2][en - 3][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            em.internal_1x2[en2b][en1b][enb][stb][st1b][st2b][r[st + 3]] +
                dp[st + 3][en - 2][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            em.internal_2x2[stb][st1b][st2b][r[st + 3]][r[en - 3]][en2b][en1b][enb] +
                dp[st + 3][en - 3][DP_P]);

        // 2x3 and 3x2 loops
        const auto two_by_three = base_internal_loop + em.internal_init[5] +
            std::min(em.internal_asym, NINIO_MAX_ASYM) +
            em.internal_2x3_mismatch[stb][st1b][en1b][enb];
        mins[DP_P] = std::min(mins[DP_P],
            two_by_three + em.InternalLoopAuGuPenalty(r[st + 3], r[en - 4]) +
                em.internal_2x3_mismatch[r[en - 4]][r[en - 3]][st2b][r[st + 3]] +
                dp[st + 3][en - 4][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            two_by_three + em.InternalLoopAuGuPenalty(r[st + 4], r[en - 3]) +
                em.internal_2x3_mismatch[r[en - 3]][r[en - 2]][r[st + 3]][r[st + 4]] +
                dp[st + 4][en - 3][DP_P]);

        // For the rest of the loops we need to apply the "other" type mismatches.
        base_internal_loop += em.internal_other_mismatch[stb][st1b][en1b][enb];

        // Lyngso for the rest.
        for (int l = 6; l <= max_inter; ++l)
          mins[DP_P] = std::min(mins[DP_P],
              lyngso[st + 2][en - 2][l - 4] - em.internal_init[l - 4] + em.internal_init[l] +
                  base_internal_loop);

        // Hairpin loops.
        mins[DP_P] = std::min(mins[DP_P], pc.Hairpin(st, en));

        const auto base_branch_cost = pc.augubranch[stb][enb] + em.multiloop_hack_a;
        // (<   ><   >)
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + dp[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        mins[DP_P] = std::min(
            mins[DP_P], base_branch_cost + dp[st + 2][en - 1][DP_U2] + em.dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        mins[DP_P] = std::min(
            mins[DP_P], base_branch_cost + dp[st + 1][en - 2][DP_U2] + em.dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        mins[DP_P] = std::min(mins[DP_P],
            base_branch_cost + dp[st + 2][en - 2][DP_U2] + em.terminal[stb][st1b][en1b][enb]);

        // (.(   ).   ) Left right coax
        for (auto cand : cand_st[CAND_P_MISMATCH])
          mins[DP_P] =
              std::min(mins[DP_P], base_branch_cost + cand.energy + dp[cand.idx + 1][en - 1][DP_U]);
        // (.(   )   .) Left outer coax
        const auto outer_coax = em.MismatchCoaxial(stb, st1b, en1b, enb);
        for (auto cand : cand_st[CAND_P_OUTER])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_mismatch_coax + outer_coax +
                  dp[cand.idx + 1][en - 2][DP_U]);
        // ((   )   ) Left flush coax
        for (auto cand : cand_st[CAND_P_FLUSH])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_flush_coax +
                  em.stack[stb][st1b][r[cand.idx]][enb] + dp[cand.idx + 1][en - 1][DP_U]);

        // (   .(   ).) Right left coax
        for (auto cand : p_cand_en[CAND_EN_P_MISMATCH][en])
          mins[DP_P] =
              std::min(mins[DP_P], base_branch_cost + cand.energy + dp[st + 1][cand.idx - 1][DP_U]);
        // (.   (   ).) Right outer coax
        for (auto cand : p_cand_en[CAND_EN_P_OUTER][en])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_mismatch_coax + outer_coax +
                  dp[st + 2][cand.idx - 1][DP_U]);
        // (   (   )) Right flush coax
        for (auto cand : p_cand_en[CAND_EN_P_FLUSH][en])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - pc.min_flush_coax +
                  em.stack[stb][r[cand.idx]][en1b][enb] + dp[st + 1][cand.idx - 1][DP_U]);

        dp[st][en][DP_P] = mins[DP_P];
      }
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        mins[DP_U] = std::min(mins[DP_U], dp[st + 1][en][DP_U]);
        mins[DP_U2] = std::min(mins[DP_U2], dp[st + 1][en][DP_U2]);
      }
      for (auto cand : cand_st[CAND_U]) {
        mins[DP_U] = std::min(mins[DP_U], cand.energy + std::min(dp[cand.idx + 1][en][DP_U], 0));
        mins[DP_U2] = std::min(mins[DP_U2], cand.energy + dp[cand.idx + 1][en][DP_U]);
      }
      for (auto cand : cand_st[CAND_U_LCOAX]) {
        const auto val =
            cand.energy + std::min(dp[cand.idx + 1][en][DP_U_WC], dp[cand.idx + 1][en][DP_U_GU]);
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_RCOAX_FWD]) {
        const auto val = cand.energy - pc.min_mismatch_coax + dp[cand.idx + 1][en][DP_U_RCOAX];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_WC_FLUSH]) {
        // (   )(<   ) > Flush coax - U
        const auto val = cand.energy + dp[cand.idx + 1][en][DP_U_WC];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_GU_FLUSH]) {
        const auto val = cand.energy + dp[cand.idx + 1][en][DP_U_GU];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_WC])
        mins[DP_U_WC] =
            std::min(mins[DP_U_WC], cand.energy + std::min(dp[cand.idx + 1][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_GU])
        mins[DP_U_GU] =
            std::min(mins[DP_U_GU], cand.energy + std::min(dp[cand.idx + 1][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_RCOAX]) {
        // (   )<.( * ). > Right coax backward
        mins[DP_U_RCOAX] =
            std::min(mins[DP_U_RCOAX], cand.energy + std::min(dp[cand.idx + 1][en][DP_U], 0));
      }

      dp[st][en][DP_U] = mins[DP_U];
      dp[st][en][DP_U2] = mins[DP_U2];
      dp[st][en][DP_U_WC] = mins[DP_U_WC];
      dp[st][en][DP_U_GU] = mins[DP_U_GU];
      dp[st][en][DP_U_RCOAX] = mins[DP_U_RCOAX];

      Energy cand_st_mins[] = {
          MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(
          sizeof(cand_st_mins) / sizeof(cand_st_mins[0]) == CAND_SIZE, "array wrong size");

      // (   ) - Normal - U, U2
      const auto normal_base = dp[st][en][DP_P] + pc.augubranch[stb][enb];
      if (normal_base < dp[st][en][DP_U] && normal_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = normal_base;

      if (IsGu(stb, enb)) {
        if (normal_base < dp[st][en][DP_U_GU] && normal_base < cand_st_mins[CAND_U_GU])
          cand_st_mins[CAND_U_GU] = normal_base;
        // Base case.
        dp[st][en][DP_U_GU] = std::min(dp[st][en][DP_U_GU], normal_base);
      } else {
        if (normal_base < dp[st][en][DP_U_WC] && normal_base < cand_st_mins[CAND_U_WC])
          cand_st_mins[CAND_U_WC] = normal_base;
        // Base case.
        dp[st][en][DP_U_WC] = std::min(dp[st][en][DP_U_WC], normal_base);
      }

      // (   ). - 3' - U, U2
      const auto dangle3_base =
          dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + em.dangle3[en1b][enb][stb];
      if (dangle3_base < dp[st][en][DP_U] && dangle3_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = dangle3_base;
      // .(   ) - 5' - U, U2
      const auto dangle5_base =
          dp[st + 1][en][DP_P] + pc.augubranch[st1b][enb] + em.dangle5[enb][stb][st1b];
      if (dangle5_base < dp[st][en][DP_U] && dangle5_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = dangle5_base;
      // .(   ). - Terminal mismatch - U, U2
      const auto terminal_base =
          dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] + em.terminal[en1b][enb][stb][st1b];
      if (terminal_base < dp[st][en][DP_U] && terminal_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = terminal_base;
      // .(   ).<(   ) > - Left coax - U, U2
      const auto lcoax_base = dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] +
          em.MismatchCoaxial(en1b, enb, stb, st1b);
      if (lcoax_base < dp[st][en][DP_U]) cand_st[CAND_U_LCOAX].push_back({lcoax_base, en});
      // (   )<.(   ). > Right coax forward - U, U2
      const auto rcoaxf_base = dp[st][en][DP_P] + pc.augubranch[stb][enb] + pc.min_mismatch_coax;
      if (rcoaxf_base < dp[st][en][DP_U]) cand_st[CAND_U_RCOAX_FWD].push_back({rcoaxf_base, en});

      // (   )<.( * ). > Right coax backward - RCOAX
      const auto rcoaxb_base = dp[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] +
          em.MismatchCoaxial(en1b, enb, stb, st1b);
      if (rcoaxb_base < dp[st][en][DP_U_RCOAX] && rcoaxb_base < cand_st_mins[CAND_U_RCOAX])
        cand_st_mins[CAND_U_RCOAX] = rcoaxb_base;
      // Base case.
      dp[st][en][DP_U_RCOAX] = std::min(dp[st][en][DP_U_RCOAX], rcoaxb_base);

      // (   )(<   ) > Flush coax - U, U2
      const auto wc_flush_base =
          dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + em.stack[en1b][enb][enb ^ 3][stb];
      const auto gu_flush_base =
          dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + em.stack[en1b][enb][enb ^ 1][stb];
      if (wc_flush_base < CAP_E && wc_flush_base < dp[st][en - 1][DP_U])
        cand_st[CAND_U_WC_FLUSH].push_back({wc_flush_base, en - 1});
      if (gu_flush_base < CAP_E && (enb == G || enb == U) && gu_flush_base < dp[st][en - 1][DP_U])
        cand_st[CAND_U_GU_FLUSH].push_back({gu_flush_base, en - 1});

      // Base cases.
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], normal_base);
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], dangle3_base);
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], dangle5_base);
      dp[st][en][DP_U] = std::min(dp[st][en][DP_U], terminal_base);
      // Note we don't include the stacking here since they can't be base cases for U.

      // Paired cases
      // (.(   )   .) Left outer coax - P
      const auto plocoax_base =
          dp[st + 2][en][DP_P] + pc.augubranch[st2b][enb] + pc.min_mismatch_coax;
      if (plocoax_base < dp[st + 1][en][DP_U]) cand_st[CAND_P_OUTER].push_back({plocoax_base, en});
      // (.   (   ).) Right outer coax
      const auto procoax_base =
          dp[st][en - 2][DP_P] + pc.augubranch[stb][en2b] + pc.min_mismatch_coax;
      if (procoax_base < dp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_OUTER][en].push_back({procoax_base, st});
      // (.(   ).   ) Left right coax
      const auto plrcoax_base = dp[st + 2][en - 1][DP_P] + pc.augubranch[st2b][en1b] +
          em.MismatchCoaxial(en1b, enb, st1b, st2b);
      if (plrcoax_base < dp[st + 1][en][DP_U])
        cand_st[CAND_P_MISMATCH].push_back({plrcoax_base, en});
      // (   .(   ).) Right left coax
      const auto prlcoax_base = dp[st + 1][en - 2][DP_P] + pc.augubranch[st1b][en2b] +
          em.MismatchCoaxial(en2b, en1b, stb, st1b);
      if (prlcoax_base < dp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_MISMATCH][en].push_back({prlcoax_base, st});
      // ((   )   ) Left flush coax
      const auto plfcoax_base = dp[st + 1][en][DP_P] + pc.augubranch[st1b][enb] + pc.min_flush_coax;
      if (plfcoax_base < dp[st + 1][en][DP_U]) cand_st[CAND_P_FLUSH].push_back({plfcoax_base, en});
      // (   (   )) Right flush coax
      const auto prfcoax_base = dp[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + pc.min_flush_coax;
      if (prfcoax_base < dp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_FLUSH][en].push_back({prfcoax_base, st});

      // Add potentials to the candidate lists.
      for (int i = 0; i < CAND_SIZE; ++i) {
        if (cand_st_mins[i] < CAP_E &&
            (cand_st[i].empty() || cand_st_mins[i] < cand_st[i].back().energy))
          cand_st[i].push_back({cand_st_mins[i], en});
      }
    }
  }

  return dp;
}

}  // namespace mrna::mfe
