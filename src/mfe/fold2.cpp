// Copyright 2016 Eliot Courtney.
#include "energy/energy_globals.h"
#include "energy/fast_energy.h"
#include "mfe/fold.h"

namespace mrna::fold::internal {

using energy::FastHairpin;
using energy::FastTwoLoop;
using energy::gem;
using energy::gpc;
using energy::ViableFoldingPair;

void ComputeTables2() {
  const int N = static_cast<int>(gr.size());
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  std::vector<std::vector<cand_t>> p_cand_en[CAND_EN_SIZE];
  for (auto& i : p_cand_en) i.resize(gr.size());
  std::vector<cand_t> cand_st[CAND_SIZE];
  for (int st = N - 1; st >= 0; --st) {
    for (auto& i : cand_st) i.clear();
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const base_t stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en],
                   en1b = gr[en - 1], en2b = gr[en - 2];
      energy_t mins[] = {MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(sizeof(mins) / sizeof(mins[0]) == DP_SIZE, "array wrong size");

      // Update paired - only if can actually pair.
      if (ViableFoldingPair(st, en)) {
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        mins[DP_P] =
            std::min(mins[DP_P], gem.stack[stb][st1b][en1b][enb] + gdp[st + 1][en - 1][DP_P]);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (gdp[ist][ien][DP_P] < mins[DP_P] - gpc.min_twoloop_not_stack)
              mins[DP_P] =
                  std::min(mins[DP_P], FastTwoLoop(st, en, ist, ien) + gdp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        mins[DP_P] = std::min(mins[DP_P], FastHairpin(st, en));

        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
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
          mins[DP_P] =
              std::min(mins[DP_P], base_branch_cost + cand.energy + gdp[cand.idx][en - 1][DP_U]);
        // (.(   )   .) Left outer coax
        const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);
        for (auto cand : cand_st[CAND_P_OUTER])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - gpc.min_mismatch_coax + outer_coax +
                  gdp[cand.idx][en - 2][DP_U]);
        // ((   )   ) Left flush coax
        for (auto cand : cand_st[CAND_P_FLUSH])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - gpc.min_flush_coax +
                  gem.stack[stb][st1b][gr[cand.idx]][enb] + gdp[cand.idx + 1][en - 1][DP_U]);
        // (   .(   ).) Right left coax
        for (auto cand : p_cand_en[CAND_EN_P_MISMATCH][en])
          mins[DP_P] =
              std::min(mins[DP_P], base_branch_cost + cand.energy + gdp[st + 1][cand.idx][DP_U]);
        // (.   (   ).) Right outer coax
        for (auto cand : p_cand_en[CAND_EN_P_OUTER][en])
          mins[DP_P] = std::min(mins[DP_P],
              base_branch_cost + cand.energy - gpc.min_mismatch_coax + outer_coax +
                  gdp[st + 2][cand.idx][DP_U]);
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
        mins[DP_U] = std::min(mins[DP_U], cand.energy + std::min(gdp[cand.idx][en][DP_U], 0));
        mins[DP_U2] = std::min(mins[DP_U2], cand.energy + gdp[cand.idx][en][DP_U]);
      }
      for (auto cand : cand_st[CAND_U_LCOAX]) {
        const auto val =
            cand.energy + std::min(gdp[cand.idx][en][DP_U_WC], gdp[cand.idx][en][DP_U_GU]);
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_RCOAX_FWD]) {
        const auto val = cand.energy - gpc.min_mismatch_coax + gdp[cand.idx][en][DP_U_RCOAX];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_WC_FLUSH]) {
        // (   )(<   ) > Flush coax - U
        const auto val = cand.energy + gdp[cand.idx][en][DP_U_WC];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_GU_FLUSH]) {
        auto val = cand.energy + gdp[cand.idx][en][DP_U_GU];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_WC])
        mins[DP_U_WC] = std::min(mins[DP_U_WC], cand.energy + std::min(gdp[cand.idx][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_GU])
        mins[DP_U_GU] = std::min(mins[DP_U_GU], cand.energy + std::min(gdp[cand.idx][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_RCOAX]) {
        // (   )<.( * ). > Right coax backward
        mins[DP_U_RCOAX] =
            std::min(mins[DP_U_RCOAX], cand.energy + std::min(gdp[cand.idx][en][DP_U], 0));
      }

      // Set these so we can use sparse folding.
      gdp[st][en][DP_U] = mins[DP_U];
      gdp[st][en][DP_U2] = mins[DP_U2];
      gdp[st][en][DP_U_WC] = mins[DP_U_WC];
      gdp[st][en][DP_U_GU] = mins[DP_U_GU];
      gdp[st][en][DP_U_RCOAX] = mins[DP_U_RCOAX];

      // Now build the candidates arrays based off the current area [st, en].
      // In general, the idea is to see if there exists something we could replace a structure with
      // that is as good.
      // e.g. we could replace a (   )3' with the equivalent U since we know the energy for (   )3'
      // that is, it is self contained. If replacing it with U[st][en] is better, then we do not
      // need to consider (...)3' when computing a larger U. In some cases we use the minimum
      // possible energy if we don't know the energy exactly for a structure (e.g. RCOAX). These
      // orderings are useful to remember: U <= U_WC, U_GU, U2
      energy_t cand_st_u = MAX_E;

      // Unpaired cases. These store the best pairs u_cand that begin at st.
      // begin means that the whole interaction starts at st. e.g. .(   ). starts one before the
      // paren.
      // (   ) - Normal - U, U2
      const auto normal_base = gdp[st][en][DP_P] + gpc.augubranch[stb][enb];
      if (normal_base < gdp[st][en][DP_U] && normal_base < cand_st_u) cand_st_u = normal_base;

      // For U_GU and U_WC, they can't be replaced with DP_U, so we need to compare them to
      // something they can be
      // replaced with, i.e. themselves.
      if (IsGu(stb, enb)) {
        if (normal_base < CAP_E && normal_base < gdp[st][en][DP_U_GU] &&
            (cand_st[CAND_U_GU].empty() || normal_base < cand_st[CAND_U_GU].back().energy))
          cand_st[CAND_U_GU].push_back({normal_base, en + 1});
        // Base case.
        gdp[st][en][DP_U_GU] = std::min(gdp[st][en][DP_U_GU], normal_base);
      } else {
        if (normal_base < CAP_E && normal_base < gdp[st][en][DP_U_WC] &&
            (cand_st[CAND_U_WC].empty() || normal_base < cand_st[CAND_U_WC].back().energy))
          cand_st[CAND_U_WC].push_back({normal_base, en + 1});
        // Base case.
        gdp[st][en][DP_U_WC] = std::min(gdp[st][en][DP_U_WC], normal_base);
      }

      // Can only merge candidate lists for monotonicity if
      // the right part of the pivot is the same (from the same array).
      // Can only apply monotonicity optimisation to ones ending with min(U, 0).
      // (   ). - 3' - U, U2
      const auto dangle3_base =
          gdp[st][en - 1][DP_P] + gpc.augubranch[stb][en1b] + gem.dangle3[en1b][enb][stb];
      if (dangle3_base < gdp[st][en][DP_U] && dangle3_base < cand_st_u) cand_st_u = dangle3_base;
      // .(   ) - 5' - U, U2
      const auto dangle5_base =
          gdp[st + 1][en][DP_P] + gpc.augubranch[st1b][enb] + gem.dangle5[enb][stb][st1b];
      if (dangle5_base < gdp[st][en][DP_U] && dangle5_base < cand_st_u) cand_st_u = dangle5_base;
      // .(   ). - Terminal mismatch - U, U2
      const auto terminal_base = gdp[st + 1][en - 1][DP_P] + gpc.augubranch[st1b][en1b] +
          gem.terminal[en1b][enb][stb][st1b];
      if (terminal_base < gdp[st][en][DP_U] && terminal_base < cand_st_u) cand_st_u = terminal_base;

      // Add potentials to the candidate lists.
      if (cand_st_u < CAP_E &&
          (cand_st[CAND_U].empty() || cand_st_u < cand_st[CAND_U].back().energy))
        cand_st[CAND_U].push_back({cand_st_u, en + 1});

      // .(   ).<(   ) > - Left coax - U, U2
      const auto lcoax_base = gdp[st + 1][en - 1][DP_P] + gpc.augubranch[st1b][en1b] +
          gem.MismatchCoaxial(en1b, enb, stb, st1b);
      if (lcoax_base < CAP_E && lcoax_base < gdp[st][en][DP_U])
        cand_st[CAND_U_LCOAX].push_back({lcoax_base, en + 1});
      // (   )<.(   ). > Right coax forward - U, U2
      const auto rcoaxf_base = gdp[st][en][DP_P] + gpc.augubranch[stb][enb] + gpc.min_mismatch_coax;
      if (rcoaxf_base < CAP_E && rcoaxf_base < gdp[st][en][DP_U])
        cand_st[CAND_U_RCOAX_FWD].push_back({rcoaxf_base, en + 1});

      // (   )<.( * ). > Right coax backward - RCOAX
      // Again, we can't replace RCOAX with U, we'd have to replace it with RCOAX, so compare to
      // itself.
      const auto rcoaxb_base = gdp[st + 1][en - 1][DP_P] + gpc.augubranch[st1b][en1b] +
          gem.MismatchCoaxial(en1b, enb, stb, st1b);
      if (rcoaxb_base < CAP_E && rcoaxb_base < gdp[st][en][DP_U_RCOAX] &&
          (cand_st[CAND_U_RCOAX].empty() || rcoaxb_base < cand_st[CAND_U_RCOAX].back().energy))
        cand_st[CAND_U_RCOAX].push_back({rcoaxb_base, en + 1});
      // Base case.
      gdp[st][en][DP_U_RCOAX] = std::min(gdp[st][en][DP_U_RCOAX], rcoaxb_base);

      // (   )(<   ) > Flush coax - U, U2
      const auto wc_flush_base =
          gdp[st][en - 1][DP_P] + gpc.augubranch[stb][en1b] + gem.stack[en1b][enb][enb ^ 3][stb];
      const auto gu_flush_base =
          gdp[st][en - 1][DP_P] + gpc.augubranch[stb][en1b] + gem.stack[en1b][enb][enb ^ 1][stb];
      if (wc_flush_base < CAP_E && wc_flush_base < gdp[st][en - 1][DP_U])
        cand_st[CAND_U_WC_FLUSH].push_back({wc_flush_base, en});
      if (gu_flush_base < CAP_E && (enb == G || enb == U) && gu_flush_base < gdp[st][en - 1][DP_U])
        cand_st[CAND_U_GU_FLUSH].push_back({gu_flush_base, en});

      // Base cases.
      gdp[st][en][DP_U] = std::min(gdp[st][en][DP_U], normal_base);
      gdp[st][en][DP_U] = std::min(gdp[st][en][DP_U], dangle3_base);
      gdp[st][en][DP_U] = std::min(gdp[st][en][DP_U], dangle5_base);
      gdp[st][en][DP_U] = std::min(gdp[st][en][DP_U], terminal_base);
      // Note we don't include the stacking here since they can't be base cases for U.

      // Paired cases
      // (.(   )   .) Left outer coax - P
      // Since we assumed the minimum energy coax stack and made this structure self contained,
      // we could potentially replace it with U[st + 1][en].
      const auto plocoax_base =
          gdp[st + 2][en][DP_P] + gpc.augubranch[st2b][enb] + gpc.min_mismatch_coax;
      if (plocoax_base < CAP_E && plocoax_base < gdp[st + 1][en][DP_U])
        cand_st[CAND_P_OUTER].push_back({plocoax_base, en + 1});
      // (.   (   ).) Right outer coax
      const auto procoax_base =
          gdp[st][en - 2][DP_P] + gpc.augubranch[stb][en2b] + gpc.min_mismatch_coax;
      if (procoax_base < CAP_E && procoax_base < gdp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_OUTER][en].push_back({procoax_base, st - 1});
      // (.(   ).   ) Left right coax
      const auto plrcoax_base = gdp[st + 2][en - 1][DP_P] + gpc.augubranch[st2b][en1b] +
          gem.MismatchCoaxial(en1b, enb, st1b, st2b);
      if (plrcoax_base < CAP_E && plrcoax_base < gdp[st + 1][en][DP_U])
        cand_st[CAND_P_MISMATCH].push_back({plrcoax_base, en + 1});
      // (   .(   ).) Right left coax
      const auto prlcoax_base = gdp[st + 1][en - 2][DP_P] + gpc.augubranch[st1b][en2b] +
          gem.MismatchCoaxial(en2b, en1b, stb, st1b);
      if (prlcoax_base < CAP_E && prlcoax_base < gdp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_MISMATCH][en].push_back({prlcoax_base, st - 1});
      // ((   )   ) Left flush coax
      const auto plfcoax_base =
          gdp[st + 1][en][DP_P] + gpc.augubranch[st1b][enb] + gpc.min_flush_coax;
      if (plfcoax_base < CAP_E && plfcoax_base < gdp[st + 1][en][DP_U])
        cand_st[CAND_P_FLUSH].push_back({plfcoax_base, en});
      // (   (   )) Right flush coax
      const auto prfcoax_base =
          gdp[st][en - 1][DP_P] + gpc.augubranch[stb][en1b] + gpc.min_flush_coax;
      if (prfcoax_base < CAP_E && prfcoax_base < gdp[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_FLUSH][en].push_back({prfcoax_base, st});
    }
  }
}

}  // namespace mrna::fold::internal
