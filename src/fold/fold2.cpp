#include "fold/fold.h"

namespace memerna {
namespace fold {

using namespace constants;
using namespace energy;
using namespace internal;

void Context::ComputeTables2() {
  static_assert(HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");

  std::vector<std::vector<cand_t>> p_cand_en[CAND_EN_SIZE];
  for (auto& i : p_cand_en) i.resize(r.size());
  std::vector<cand_t> cand_st[CAND_SIZE];
  for (int st = N - 1; st >= 0; --st) {
    for (int i = 0; i < CAND_SIZE; ++i)
      cand_st[i].clear();
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      base_t stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];
      energy_t mins[] = {MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(sizeof(mins) / sizeof(mins[0]) == DP_SIZE, "array wrong size");

      // Update paired - only if can actually pair.
      if (ViableFoldingPair(r, st, en)) {
        int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        mins[DP_P] = std::min(mins[DP_P], em.stack[stb][st1b][en1b][enb] + arr[st + 1][en - 1][DP_P]);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (arr[ist][ien][DP_P] < mins[DP_P] - pc.min_twoloop_not_stack)
              mins[DP_P] = std::min(mins[DP_P], FastTwoLoop(st, en, ist, ien) + arr[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        mins[DP_P] = std::min(mins[DP_P], FastHairpin(st, en));

        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = pc.augubranch[stb][enb] + em.multiloop_hack_a;

        // (<   ><   >)
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + arr[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + arr[st + 2][en - 1][DP_U2] + em.dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + arr[st + 1][en - 2][DP_U2] + em.dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        mins[DP_P] = std::min(mins[DP_P],
            base_branch_cost + arr[st + 2][en - 2][DP_U2] + em.terminal[stb][st1b][en1b][enb]);

        // (.(   ).   ) Left right coax
        for (auto cand : cand_st[CAND_P_MISMATCH])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy + arr[cand.idx + 1][en - 1][DP_U]);
        // (.(   )   .) Left outer coax
        auto outer_coax = em.MismatchCoaxial(stb, st1b, en1b, enb);
        for (auto cand : cand_st[CAND_P_OUTER])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy - pc.min_mismatch_coax +
              outer_coax + arr[cand.idx + 1][en - 2][DP_U]);
        // ((   )   ) Left flush coax
        for (auto cand : cand_st[CAND_P_FLUSH])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy - pc.min_flush_coax +
              em.stack[stb][st1b][r[cand.idx]][enb] + arr[cand.idx + 1][en - 1][DP_U]);
        // (   .(   ).) Right left coax
        for (auto cand : p_cand_en[CAND_EN_P_MISMATCH][en])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy + arr[st + 1][cand.idx - 1][DP_U]);
        // (.   (   ).) Right outer coax
        for (auto cand : p_cand_en[CAND_EN_P_OUTER][en])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy - pc.min_mismatch_coax +
              outer_coax + arr[st + 2][cand.idx - 1][DP_U]);
        // (   (   )) Right flush coax
        for (auto cand : p_cand_en[CAND_EN_P_FLUSH][en])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy - pc.min_flush_coax +
              em.stack[stb][r[cand.idx]][en1b][enb] + arr[st + 1][cand.idx - 1][DP_U]);

        arr[st][en][DP_P] = mins[DP_P];
      }
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        mins[DP_U] = std::min(mins[DP_U], arr[st + 1][en][DP_U]);
        mins[DP_U2] = std::min(mins[DP_U2], arr[st + 1][en][DP_U2]);
      }
      for (auto cand : cand_st[CAND_U]) {
        mins[DP_U] = std::min(mins[DP_U], cand.energy + std::min(arr[cand.idx + 1][en][DP_U], 0));
        mins[DP_U2] = std::min(mins[DP_U2], cand.energy + arr[cand.idx + 1][en][DP_U]);
      }
      for (auto cand : cand_st[CAND_U_LCOAX]) {
        auto val = cand.energy + std::min(arr[cand.idx + 1][en][DP_U_WC], arr[cand.idx + 1][en][DP_U_GU]);
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_RCOAX_FWD]) {
        auto val = cand.energy - pc.min_mismatch_coax + arr[cand.idx + 1][en][DP_U_RCOAX];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_WC_FLUSH]) {
        // (   )<(   ) > Flush coax - U
        auto val = cand.energy + arr[cand.idx + 1][en][DP_U_WC];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_GU_FLUSH]) {
        auto val = cand.energy + arr[cand.idx + 1][en][DP_U_GU];
        mins[DP_U] = std::min(mins[DP_U], val);
        mins[DP_U2] = std::min(mins[DP_U2], val);
      }
      for (auto cand : cand_st[CAND_U_WC])
        mins[DP_U_WC] = std::min(mins[DP_U_WC], cand.energy + std::min(arr[cand.idx + 1][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_GU])
        mins[DP_U_GU] = std::min(mins[DP_U_GU], cand.energy + std::min(arr[cand.idx + 1][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_RCOAX]) {
        // (   ).<( * ). > Right coax backward
        assert(st > 0);
        mins[DP_U_RCOAX] = std::min(mins[DP_U_RCOAX], cand.energy + std::min(arr[cand.idx + 1][en][DP_U], 0));
      }

      // Set these so we can use sparse folding.
      arr[st][en][DP_U] = mins[DP_U];
      arr[st][en][DP_U2] = mins[DP_U2];
      arr[st][en][DP_U_WC] = mins[DP_U_WC];
      arr[st][en][DP_U_GU] = mins[DP_U_GU];
      arr[st][en][DP_U_RCOAX] = mins[DP_U_RCOAX];

      // Now build the candidates arrays based off the current area [st, en].
      // In general, the idea is to see if there exists something we could replace a structure with that is as good.
      // e.g. we could replace a (   )3' with the equivalent U since we know the energy for (   )3' -
      // it is self contained. If replacing it with U[st][en] is better, then we do not need to consider (...)3'
      // when computing a larger U. In some cases we use the minimum possible energy if we don't know the energy exactly
      // for a structure (e.g. RCOAX).
      // These orderings are useful to remember:
      // U <= U_WC, U_GU, U2

      energy_t cand_st_mins[] = {MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(sizeof(cand_st_mins) / sizeof(cand_st_mins[0]) == CAND_SIZE, "array wrong size");

      // Unpaired cases. These store the best pairs u_cand that begin at st.
      // begin means that the whole interaction starts at st. e.g. .(   ). starts one before the paren.
      // (   ) - Normal - U, U2
      auto normal_base = arr[st][en][DP_P] + pc.augubranch[stb][enb];
      if (normal_base < arr[st][en][DP_U] && normal_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = normal_base;

      // For U_GU and U_WC, they can't be replaced with DP_U, so we need to compare them to something they can be
      // replaced with, i.e. themselves.
      if (IsGu(stb, enb)) {
        if (normal_base < arr[st][en][DP_U_GU] && normal_base < cand_st_mins[CAND_U_GU])
          cand_st_mins[CAND_U_GU] = normal_base;
        // Base case.
        arr[st][en][DP_U_GU] = std::min(arr[st][en][DP_U_GU], normal_base);
      } else {
        if (normal_base < arr[st][en][DP_U_WC] && normal_base < cand_st_mins[CAND_U_WC])
          cand_st_mins[CAND_U_WC] = normal_base;
        // Base case.
        arr[st][en][DP_U_WC] = std::min(arr[st][en][DP_U_WC], normal_base);
      }

      // Can only merge candidate lists for monotonicity if
      // the right part of the pivot is the same (from the same array).
      // Can only apply monotonicity optimisation to ones ending with min(U, 0).
      // (   ). - 3' - U, U2
      auto dangle3_base = arr[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + em.dangle3[en1b][enb][stb];
      if (dangle3_base < arr[st][en][DP_U] && dangle3_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = dangle3_base;
      // .(   ) - 5' - U, U2
      auto dangle5_base = arr[st + 1][en][DP_P] + pc.augubranch[st1b][enb] + em.dangle5[enb][stb][st1b];
      if (dangle5_base < arr[st][en][DP_U] && dangle5_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = dangle5_base;
      // .(   ). - Terminal mismatch - U, U2
      auto terminal_base = arr[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] + em.terminal[en1b][enb][stb][st1b];
      if (terminal_base < arr[st][en][DP_U] && terminal_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = terminal_base;
      // .(   ).<(   ) > - Left coax - U, U2
      auto lcoax_base = arr[st + 1][en - 1][DP_P] + pc.augubranch[st1b][en1b] +
          em.MismatchCoaxial(en1b, enb, stb, st1b);
      if (lcoax_base < CAP_E && lcoax_base < arr[st][en][DP_U])
        cand_st[CAND_U_LCOAX].push_back({lcoax_base, en});
      // (   ).<(   ). > Right coax forward - U, U2
      // This is probably better than having four candidate lists for each possible mismatch (TODO check this).
      auto rcoaxf_base = arr[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + pc.min_mismatch_coax;
      if (rcoaxf_base < CAP_E && rcoaxf_base < arr[st][en][DP_U])
        cand_st[CAND_U_RCOAX_FWD].push_back({rcoaxf_base, en});

      // (   ).<( * ). > Right coax backward - RCOAX
      // Again, we can't replace RCOAX with U, we'd have to replace it with RCOAX, so compare to itself.
      if (st > 0) {
        auto rcoaxb_base = arr[st][en - 1][DP_P] + pc.augubranch[stb][en1b] +
            em.MismatchCoaxial(en1b, enb, r[st - 1], stb);
        if (rcoaxb_base < arr[st][en][DP_U_RCOAX] && rcoaxb_base < cand_st_mins[CAND_U_RCOAX])
          cand_st_mins[CAND_U_RCOAX] = rcoaxb_base;
        // Base case.
        arr[st][en][DP_U_RCOAX] = std::min(arr[st][en][DP_U_RCOAX], rcoaxb_base);
      }

      // (   )<(   ) > Flush coax - U, U2
      if (en + 1 < N) {
        auto enr1b = r[en + 1];
        auto wc_flush_base = arr[st][en][DP_P] + pc.augubranch[stb][enb] + em.stack[enb][enr1b][enr1b ^ 3][stb];
        auto gu_flush_base = arr[st][en][DP_P] + pc.augubranch[stb][enb] + em.stack[enb][enr1b][enr1b ^ 1][stb];
        if (wc_flush_base < CAP_E && wc_flush_base < arr[st][en][DP_U])
          cand_st[CAND_U_WC_FLUSH].push_back({wc_flush_base, en});
        if (gu_flush_base < CAP_E && (enr1b == G || enr1b == U) && gu_flush_base < arr[st][en][DP_U])
          cand_st[CAND_U_GU_FLUSH].push_back({gu_flush_base, en});
      }

      // Base cases.
      arr[st][en][DP_U] = std::min(arr[st][en][DP_U], normal_base);
      arr[st][en][DP_U] = std::min(arr[st][en][DP_U], dangle3_base);
      arr[st][en][DP_U] = std::min(arr[st][en][DP_U], dangle5_base);
      arr[st][en][DP_U] = std::min(arr[st][en][DP_U], terminal_base);
      // Note we don't include the stacking here since they can't be base cases for U.

      // Paired cases
      // TODO is it better to use more arrays for each possible coax stack, or min coax here?
      // TODO could also use the one pair we know to find out the min poss stack -- test if perf gain useful
      // (.(   )   .) Left outer coax - P
      // Since we assumed the minimum energy coax stack and made this structure self contained,
      // we could potentially replace it with U[st + 1][en].
      auto plocoax_base = arr[st + 2][en][DP_P] + pc.augubranch[st2b][enb] + pc.min_mismatch_coax;
      if (plocoax_base < CAP_E && plocoax_base < arr[st + 1][en][DP_U])
        cand_st[CAND_P_OUTER].push_back({plocoax_base, en});
      // (.   (   ).) Right outer coax
      auto procoax_base = arr[st][en - 2][DP_P] + pc.augubranch[stb][en2b] + pc.min_mismatch_coax;
      if (procoax_base < CAP_E && procoax_base < arr[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_OUTER][en].push_back({procoax_base, st});
      // (.(   ).   ) Left right coax
      auto plrcoax_base = arr[st + 2][en - 1][DP_P] + pc.augubranch[st2b][en1b] +
          em.MismatchCoaxial(en1b, enb, st1b, st2b);
      if (plrcoax_base < CAP_E && plrcoax_base < arr[st + 1][en][DP_U])
        cand_st[CAND_P_MISMATCH].push_back({plrcoax_base, en});
      // (   .(   ).) Right left coax
      auto prlcoax_base = arr[st + 1][en - 2][DP_P] + pc.augubranch[st1b][en2b] +
          em.MismatchCoaxial(en2b, en1b, stb, st1b);
      if (prlcoax_base < CAP_E && prlcoax_base < arr[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_MISMATCH][en].push_back({prlcoax_base, st});
      // ((   )   ) Left flush coax
      auto plfcoax_base = arr[st + 1][en][DP_P] + pc.augubranch[st1b][enb] + pc.min_flush_coax;
      if (plfcoax_base < CAP_E && plfcoax_base < arr[st + 1][en][DP_U])
        cand_st[CAND_P_FLUSH].push_back({plfcoax_base, en});
      // (   (   )) Right flush coax
      auto prfcoax_base = arr[st][en - 1][DP_P] + pc.augubranch[stb][en1b] + pc.min_flush_coax;
      if (prfcoax_base < CAP_E && prfcoax_base < arr[st][en - 1][DP_U])
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

}
}
