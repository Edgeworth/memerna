// Copyright 2016 Eliot Courtney.
#include "energy/energy_globals.h"
#include "energy/energy_model.h"
#include "fold/fold.h"

namespace mrna {
namespace fold {
namespace internal {

using energy::EnergyModel;
using energy::FastHairpin;
using energy::gem;
using energy::gpc;
using energy::ViableFoldingPair;

void ComputeTables3() {
  const int N = static_cast<int>(gr.size());
  static_assert(
      HAIRPIN_MIN_SZ >= 3, "Minimum hairpin size >= 3 is relied upon in some expressions.");

  // See ComputeTables2 for comments - it is mostly the same.
  std::vector<std::vector<cand_t>> p_cand_en[CAND_EN_SIZE];
  for (auto& i : p_cand_en) i.resize(gr.size());
  std::vector<cand_t> cand_st[CAND_SIZE];
  array3d_t<energy_t, TWOLOOP_MAX_SZ + 1> lyngso(gr.size());
  for (int st = N - 1; st >= 0; --st) {
    for (auto& i : cand_st) i.clear();
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const base_t stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en],
                   en1b = gr[en - 1], en2b = gr[en - 2];
      energy_t mins[] = {MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(sizeof(mins) / sizeof(mins[0]) == DP_SIZE, "array wrong size");
      const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);

      for (int l = 0; l <= max_inter; ++l) {
        // Don't add asymmetry here
        if (l >= 2)
          lyngso[st][en][l] = std::min(lyngso[st][en][l],
              lyngso[st + 1][en - 1][l - 2] - gem.internal_init[l - 2] + gem.internal_init[l]);

        // Add asymmetry here, on left and right
        auto val = std::min(l * gem.internal_asym, NINIO_MAX_ASYM) + gem.internal_init[l];
        lyngso[st][en][l] = std::min(lyngso[st][en][l],
            gem.InternalLoopAuGuPenalty(gr[st + l + 1], en1b) +
                gem.internal_other_mismatch[en1b][enb][gr[st + l]][gr[st + l + 1]] + val +
                gdp[st + l + 1][en - 1][DP_P]);
        lyngso[st][en][l] = std::min(lyngso[st][en][l],
            gem.InternalLoopAuGuPenalty(st1b, gr[en - l - 1]) +
                gem.internal_other_mismatch[gr[en - l - 1]][gr[en - l]][stb][st1b] + val +
                gdp[st + 1][en - l - 1][DP_P]);
      }

      // Update paired - only if can actually pair.
      if (ViableFoldingPair(st, en)) {
        // Stacking
        mins[DP_P] =
            std::min(mins[DP_P], gem.stack[stb][st1b][en1b][enb] + gdp[st + 1][en - 1][DP_P]);
        // Bulge
        for (int isz = 1; isz <= max_inter; ++isz) {
          mins[DP_P] = std::min(mins[DP_P],
              gem.Bulge(gr, st, en, st + 1 + isz, en - 1) + gdp[st + 1 + isz][en - 1][DP_P]);
          mins[DP_P] = std::min(mins[DP_P],
              gem.Bulge(gr, st, en, st + 1, en - 1 - isz) + gdp[st + 1][en - 1 - isz][DP_P]);
        }

        // Ax1 internal loops. Make sure to skip 0x1, 1x1, 2x1, and 1x2 loops, since they have
        // special energies.
        static_assert(EnergyModel::INITIATION_CACHE_SZ > TWOLOOP_MAX_SZ,
            "need initiation cached up to TWOLOOP_MAX_SZ");
        auto base_internal_loop = gem.InternalLoopAuGuPenalty(stb, enb);
        for (int isz = 4; isz <= max_inter; ++isz) {
          auto val = base_internal_loop + gem.internal_init[isz] +
              std::min((isz - 2) * gem.internal_asym, NINIO_MAX_ASYM);
          mins[DP_P] = std::min(mins[DP_P],
              val + gem.InternalLoopAuGuPenalty(gr[st + isz], en2b) + gdp[st + isz][en - 2][DP_P]);
          mins[DP_P] = std::min(mins[DP_P],
              val + gem.InternalLoopAuGuPenalty(st2b, gr[en - isz]) + gdp[st + 2][en - isz][DP_P]);
        }

        // Internal loop cases. Since we require HAIRPIN_MIN_SZ >= 3 and initialise arr to MAX_E, we
        // don't need ifs
        // here.
        mins[DP_P] = std::min(mins[DP_P],
            gem.internal_1x1[stb][st1b][st2b][en2b][en1b][enb] + gdp[st + 2][en - 2][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            gem.internal_1x2[stb][st1b][st2b][gr[en - 3]][en2b][en1b][enb] +
                gdp[st + 2][en - 3][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            gem.internal_1x2[en2b][en1b][enb][stb][st1b][st2b][gr[st + 3]] +
                gdp[st + 3][en - 2][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            gem.internal_2x2[stb][st1b][st2b][gr[st + 3]][gr[en - 3]][en2b][en1b][enb] +
                gdp[st + 3][en - 3][DP_P]);

        // 2x3 and 3x2 loops
        const auto two_by_three = base_internal_loop + gem.internal_init[5] +
            std::min(gem.internal_asym, NINIO_MAX_ASYM) +
            gem.internal_2x3_mismatch[stb][st1b][en1b][enb];
        mins[DP_P] = std::min(mins[DP_P],
            two_by_three + gem.InternalLoopAuGuPenalty(gr[st + 3], gr[en - 4]) +
                gem.internal_2x3_mismatch[gr[en - 4]][gr[en - 3]][st2b][gr[st + 3]] +
                gdp[st + 3][en - 4][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            two_by_three + gem.InternalLoopAuGuPenalty(gr[st + 4], gr[en - 3]) +
                gem.internal_2x3_mismatch[gr[en - 3]][gr[en - 2]][gr[st + 3]][gr[st + 4]] +
                gdp[st + 4][en - 3][DP_P]);

        // For the rest of the loops we need to apply the "other" type mismatches.
        base_internal_loop += gem.internal_other_mismatch[stb][st1b][en1b][enb];

        // Lyngso for the rest.
        for (int l = 6; l <= max_inter; ++l)
          mins[DP_P] = std::min(mins[DP_P],
              lyngso[st + 2][en - 2][l - 4] - gem.internal_init[l - 4] + gem.internal_init[l] +
                  base_internal_loop);

        // Hairpin loops.
        mins[DP_P] = std::min(mins[DP_P], FastHairpin(st, en));

        const auto base_branch_cost = gpc.augubranch[stb][enb] + gem.multiloop_hack_a;
        // (<   ><   >)
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + gdp[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        mins[DP_P] = std::min(mins[DP_P],
            base_branch_cost + gdp[st + 2][en - 1][DP_U2] + gem.dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        mins[DP_P] = std::min(mins[DP_P],
            base_branch_cost + gdp[st + 1][en - 2][DP_U2] + gem.dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        mins[DP_P] = std::min(mins[DP_P],
            base_branch_cost + gdp[st + 2][en - 2][DP_U2] + gem.terminal[stb][st1b][en1b][enb]);

        // (.(   ).   ) Left right coax
        for (auto cand : cand_st[CAND_P_MISMATCH])
          mins[DP_P] = std::min(
              mins[DP_P], base_branch_cost + cand.energy + gdp[cand.idx + 1][en - 1][DP_U]);
        // (.(   )   .) Left outer coax
        const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);
        for (auto cand : cand_st[CAND_P_OUTER])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - gpc.min_mismatch_coax + outer_coax +
                  gdp[cand.idx + 1][en - 2][DP_U]);
        // ((   )   ) Left flush coax
        for (auto cand : cand_st[CAND_P_FLUSH])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - gpc.min_flush_coax +
                  gem.stack[stb][st1b][gr[cand.idx]][enb] + gdp[cand.idx + 1][en - 1][DP_U]);

        // (   .(   ).) Right left coax
        for (auto cand : p_cand_en[CAND_EN_P_MISMATCH][en])
          mins[DP_P] = std::min(
              mins[DP_P], base_branch_cost + cand.energy + gdp[st + 1][cand.idx - 1][DP_U]);
        // (.   (   ).) Right outer coax
        for (auto cand : p_cand_en[CAND_EN_P_OUTER][en])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - gpc.min_mismatch_coax + outer_coax +
                  gdp[st + 2][cand.idx - 1][DP_U]);
        // (   (   )) Right flush coax
        for (auto cand : p_cand_en[CAND_EN_P_FLUSH][en])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - gpc.min_flush_coax +
                  gem.stack[stb][gr[cand.idx]][en1b][enb] + gdp[st + 1][cand.idx - 1][DP_U]);

        gdp[st][en][DP_P] = mins[DP_P];
      }
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        mins[DP_U] = std::min(mins[DP_U], gdp[st + 1][en][DP_U]);
        mins[DP_U2] = std::min(mins[DP_U2], gdp[st + 1][en][DP_U2]);
      }
      for (auto cand : cand_st[CAND_U]) {
        mins[DP_U] = std::min(mins[DP_U], cand.energy + std::min(gdp[cand.idx + 1][en][DP_U], 0));
        mins[DP_U2] = std::min(mins[DP_U2], cand.energy + gdp[cand.idx + 1][en][DP_U]);
      }
      for (auto cand : cand_st[CAND_U_LCOAX]) {
        const auto val =
            cand.energy + std::min(gdp[cand.idx + 1][en][DP_U_WC], gdp[cand.idx + 1][en][DP_U_GU]);
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_RCOAX_FWD]) {
        const auto val = cand.energy - gpc.min_mismatch_coax + gdp[cand.idx + 1][en][DP_U_RCOAX];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_WC_FLUSH]) {
        // (   )(<   ) > Flush coax - U
        const auto val = cand.energy + gdp[cand.idx + 1][en][DP_U_WC];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_GU_FLUSH]) {
        const auto val = cand.energy + gdp[cand.idx + 1][en][DP_U_GU];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_WC])
        mins[DP_U_WC] =
            std::min(mins[DP_U_WC], cand.energy + std::min(gdp[cand.idx + 1][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_GU])
        mins[DP_U_GU] =
            std::min(mins[DP_U_GU], cand.energy + std::min(gdp[cand.idx + 1][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_RCOAX]) {
        // (   )<.( * ). > Right coax backward
        mins[DP_U_RCOAX] =
            std::min(mins[DP_U_RCOAX], cand.energy + std::min(gdp[cand.idx + 1][en][DP_U], 0));
      }

      gdp[st][en][DP_U] = mins[DP_U];
      gdp[st][en][DP_U2] = mins[DP_U2];
      gdp[st][en][DP_U_WC] = mins[DP_U_WC];
      gdp[st][en][DP_U_GU] = mins[DP_U_GU];
      gdp[st][en][DP_U_RCOAX] = mins[DP_U_RCOAX];

      energy_t cand_st_mins[] = {
          MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(
          sizeof(cand_st_mins) / sizeof(cand_st_mins[0]) == CAND_SIZE, "array wrong size");

      // (   ) - Normal - U, U2
      const auto normal_base = gdp[st][en][DP_P] + gpc.augubranch[stb][enb];
      if (normal_base < gdp[st][en][DP_U] && normal_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = normal_base;

      if (IsGu(stb, enb)) {
        if (normal_base < gdp[st][en][DP_U_GU] && normal_base < cand_st_mins[CAND_U_GU])
          cand_st_mins[CAND_U_GU] = normal_base;
        // Base case.
        gdp[st][en][DP_U_GU] = std::min(gdp[st][en][DP_U_GU], normal_base);
      } else {
        if (normal_base < gdp[st][en][DP_U_WC] && normal_base < cand_st_mins[CAND_U_WC])
          cand_st_mins[CAND_U_WC] = normal_base;
        // Base case.
        gdp[st][en][DP_U_WC] = std::min(gdp[st][en][DP_U_WC], normal_base);
      }

      // (   ). - 3' - U, U2
      const auto dangle3_base =
          gdp[st][en - 1][DP_P] + gpc.augubranch[stb][en1b] + gem.dangle3[en1b][enb][stb];
      if (dangle3_base < gdp[st][en][DP_U] && dangle3_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = dangle3_base;
      // .(   ) - 5' - U, U2
      const auto dangle5_base =
          gdp[st + 1][en][DP_P] + gpc.augubranch[st1b][enb] + gem.dangle5[enb][stb][st1b];
      if (dangle5_base < gdp[st][en][DP_U] && dangle5_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = dangle5_base;
      // .(   ). - Terminal mismatch - U, U2
      const auto terminal_base = gdp[st + 1][en - 1][DP_P] + gpc.augubranch[st1b][en1b] +
          gem.terminal[en1b][enb][stb][st1b];
      if (terminal_base < gdp[st][en][DP_U] && terminal_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = terminal_base;
      // .(   ).<(   ) > - Left coax - U, U2
      const auto lcoax_base = gdp[st + 1][en - 1][DP_P] + gpc.augubranch[st1b][en1b] +
          gem.MismatchCoaxial(en1b, enb, stb, st1b);
      if (lcoax_base < gdp[st][en][DP_U]) cand_st[CAND_U_LCOAX].push_back({lcoax_base, en});
      // (   )<.(   ). > Right coax forward - U, U2
      const auto rcoaxf_base = gdp[st][en][DP_P] + gpc.augubranch[stb][enb] + gpc.min_mismatch_coax;
      if (rcoaxf_base < gdp[st][en][DP_U]) cand_st[CAND_U_RCOAX_FWD].push_back({rcoaxf_base, en});

      // (   )<.( * ). > Right coax backward - RCOAX
      const auto rcoaxb_base = gdp[st + 1][en - 1][DP_P] + gpc.augubranch[st1b][en1b] +
          gem.MismatchCoaxial(en1b, enb, stb, st1b);
      if (rcoaxb_base < gdp[st][en][DP_U_RCOAX] && rcoaxb_base < cand_st_mins[CAND_U_RCOAX])
        cand_st_mins[CAND_U_RCOAX] = rcoaxb_base;
      // Base case.
      gdp[st][en][DP_U_RCOAX] = std::min(gdp[st][en][DP_U_RCOAX], rcoaxb_base);

      // (   )(<   ) > Flush coax - U, U2
      const auto wc_flush_base =
          gdp[st][en - 1][DP_P] + gpc.augubranch[stb][en1b] + gem.stack[en1b][enb][enb ^ 3][stb];
      const auto gu_flush_base =
          gdp[st][en - 1][DP_P] + gpc.augubranch[stb][en1b] + gem.stack[en1b][enb][enb ^ 1][stb];
      if (wc_flush_base < CAP_E && wc_flush_base < gdp[st][en - 1][DP_U])
        cand_st[CAND_U_WC_FLUSH].push_back({wc_flush_base, en - 1});
      if (gu_flush_base < CAP_E && (enb == G || enb == U) && gu_flush_base < gdp[st][en - 1][DP_U])
        cand_st[CAND_U_GU_FLUSH].push_back({gu_flush_base, en - 1});

      // Base cases.
      gdp[st][en][DP_U] = std::min(gdp[st][en][DP_U], normal_base);
      gdp[st][en][DP_U] = std::min(gdp[st][en][DP_U], dangle3_base);
      gdp[st][en][DP_U] = std::min(gdp[st][en][DP_U], dangle5_base);
      gdp[st][en][DP_U] = std::min(gdp[st][en][DP_U], terminal_base);
      // Note we don't include the stacking here since they can't be base cases for U.

      // Paired cases
      // (.(   )   .) Left outer coax - P
      const auto plocoax_base =
          gdp[st + 2][en][DP_P] + gpc.augubranch[st2b][enb] + gpc.min_mismatch_coax;
      if (plocoax_base < gdp[st + 1][en][DP_U]) cand_st[CAND_P_OUTER].push_back({plocoax_base, en});
      // (.   (   ).) Right outer coax
      const auto procoax_base =
          gdp[st][en - 2][DP_P] + gpc.augubranch[stb][en2b] + gpc.min_mismatch_coax;
      if (procoax_base < gdp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_OUTER][en].push_back({procoax_base, st});
      // (.(   ).   ) Left right coax
      const auto plrcoax_base = gdp[st + 2][en - 1][DP_P] + gpc.augubranch[st2b][en1b] +
          gem.MismatchCoaxial(en1b, enb, st1b, st2b);
      if (plrcoax_base < gdp[st + 1][en][DP_U])
        cand_st[CAND_P_MISMATCH].push_back({plrcoax_base, en});
      // (   .(   ).) Right left coax
      const auto prlcoax_base = gdp[st + 1][en - 2][DP_P] + gpc.augubranch[st1b][en2b] +
          gem.MismatchCoaxial(en2b, en1b, stb, st1b);
      if (prlcoax_base < gdp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_MISMATCH][en].push_back({prlcoax_base, st});
      // ((   )   ) Left flush coax
      const auto plfcoax_base =
          gdp[st + 1][en][DP_P] + gpc.augubranch[st1b][enb] + gpc.min_flush_coax;
      if (plfcoax_base < gdp[st + 1][en][DP_U]) cand_st[CAND_P_FLUSH].push_back({plfcoax_base, en});
      // (   (   )) Right flush coax
      const auto prfcoax_base =
          gdp[st][en - 1][DP_P] + gpc.augubranch[stb][en1b] + gpc.min_flush_coax;
      if (prfcoax_base < gdp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_FLUSH][en].push_back({prfcoax_base, st});

      // Add potentials to the candidate lists.
      for (int i = 0; i < CAND_SIZE; ++i) {
        if (cand_st_mins[i] < CAP_E &&
            (cand_st[i].empty() || cand_st_mins[i] < cand_st[i].back().energy))
          cand_st[i].push_back({cand_st_mins[i], en});
      }
    }
  }
}
}  // namespace internal
}  // namespace fold
}  // namespace mrna
