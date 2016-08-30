#include "fold/fold.h"
#include "fold/fold_globals.h"

namespace memerna {
namespace fold {

using constants::MAX_E;
using constants::CAP_E;
using constants::TWOLOOP_MAX_SZ;
using constants::NINIO_MAX_ASYM;
using constants::HAIRPIN_MIN_SZ;

array3d_t<energy_t, DP_SIZE> ComputeTables3() {
  InitFold();
  int N = int(r.size());
  // Automatically initialised to MAX_E.
  array3d_t<energy_t, DP_SIZE> arr(r.size() + 1);
  std::vector<std::vector<cand_t>> p_cand_en[CAND_EN_SIZE];
  for (auto& i : p_cand_en) i.resize(r.size());
  std::vector<cand_t> cand_st[CAND_SIZE];
  // Hairpin optimisation
  auto hairpin_precomp = PrecomputeFastHairpin();

  array3d_t<energy_t, TWOLOOP_MAX_SZ + 1> lyngso(r.size());
  static_assert(HAIRPIN_MIN_SZ >= 3, "Minimum hairpin size >= 3 is relied upon in some expressions.");
  for (int st = N - 1; st >= 0; --st) {
    for (int i = 0; i < CAND_SIZE; ++i)
      cand_st[i].clear();
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      base_t stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];
      energy_t mins[] = {MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(sizeof(mins) / sizeof(mins[0]) == DP_SIZE, "array wrong size");
      int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);

      // Lyngso for the rest.
      for (int l = 0; l <= max_inter; ++l) {
        // Don't add asymmetry here
        if (l >= 2)
          lyngso[st][en][l] = std::min(lyngso[st][en][l],
              lyngso[st + 1][en - 1][l - 2] - g_internal_init[l - 2] + g_internal_init[l]);

        // Add asymmetry here, on left and right
        auto val = std::min(l * g_internal_asym, NINIO_MAX_ASYM) + g_internal_init[l];
        lyngso[st][en][l] = std::min(lyngso[st][en][l], energy::InternalLoopAuGuPenalty(r[st + l + 1], en1b) +
            g_internal_other_mismatch[en1b][enb][r[st + l]][r[st + l + 1]] + val + arr[st + l + 1][en - 1][DP_P]);
        lyngso[st][en][l] = std::min(lyngso[st][en][l], energy::InternalLoopAuGuPenalty(st1b, r[en - l - 1]) +
            g_internal_other_mismatch[r[en - l - 1]][r[en - l]][stb][st1b] + val + arr[st + 1][en - l - 1][DP_P]);
      }

      // Update paired - only if can actually pair.
      if (CanPair(r[st], r[en]) && IsNotLonely(st, en)) {
        // Stacking
        mins[DP_P] = std::min(mins[DP_P], g_stack[stb][st1b][en1b][enb] + arr[st + 1][en - 1][DP_P]);
        // Bulge
        for (int isz = 1; isz <= max_inter; ++isz) {
          mins[DP_P] = std::min(mins[DP_P],
              energy::Bulge(st, en, st + 1 + isz, en - 1) + arr[st + 1 + isz][en - 1][DP_P]);
          mins[DP_P] = std::min(mins[DP_P],
              energy::Bulge(st, en, st + 1, en - 1 - isz) + arr[st + 1][en - 1 - isz][DP_P]);
        }

        // Ax1 internal loops. Make sure to skip 0x1, 1x1, 2x1, and 1x2 loops, since they have special energies.
        static_assert(INITIATION_CACHE_SZ > TWOLOOP_MAX_SZ, "need initiation cached up to TWOLOOP_MAX_SZ");
        auto base_internal_loop = energy::InternalLoopAuGuPenalty(stb, enb);
        for (int isz = 4; isz <= max_inter; ++isz) {
          auto val = base_internal_loop + g_internal_init[isz] + std::min((isz - 2) * g_internal_asym, NINIO_MAX_ASYM);
          mins[DP_P] = std::min(mins[DP_P], val +
              energy::InternalLoopAuGuPenalty(r[st + isz], en2b) + arr[st + isz][en - 2][DP_P]);
          mins[DP_P] = std::min(mins[DP_P], val +
              energy::InternalLoopAuGuPenalty(st2b, r[en - isz]) + arr[st + 2][en - isz][DP_P]);
        }

        // Internal loop cases. Since we require HAIRPIN_MIN_SZ >= 3 and initialise arr to MAX_E, we don't need ifs here.
        mins[DP_P] = std::min(mins[DP_P], g_internal_1x1[stb][st1b][st2b][en2b][en1b][enb] + arr[st + 2][en - 2][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            g_internal_1x2[stb][st1b][st2b][r[en - 3]][en2b][en1b][enb] + arr[st + 2][en - 3][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            g_internal_1x2[en2b][en1b][enb][stb][st1b][st2b][r[st + 3]] + arr[st + 3][en - 2][DP_P]);
        mins[DP_P] = std::min(mins[DP_P],
            g_internal_2x2[stb][st1b][st2b][r[st + 3]][r[en - 3]][en2b][en1b][enb] + arr[st + 3][en - 3][DP_P]);

        // 2x3 and 3x2 loops
        auto two_by_three = base_internal_loop + g_internal_init[5] +
            std::min(g_internal_asym, NINIO_MAX_ASYM) +
            g_internal_2x3_mismatch[stb][st1b][en1b][enb];
        mins[DP_P] = std::min(mins[DP_P], two_by_three + energy::InternalLoopAuGuPenalty(r[st + 3], r[en - 4]) +
            g_internal_2x3_mismatch[r[en - 4]][r[en - 3]][st2b][r[st + 3]] + arr[st + 3][en - 4][DP_P]);
        mins[DP_P] = std::min(mins[DP_P], two_by_three + energy::InternalLoopAuGuPenalty(r[st + 4], r[en - 3]) +
            g_internal_2x3_mismatch[r[en - 3]][r[en - 2]][r[st + 3]][r[st + 4]] + arr[st + 4][en - 3][DP_P]);

        // For the rest of the loops we need to apply the "other" type mismatches.
        base_internal_loop += g_internal_other_mismatch[stb][st1b][en1b][enb];

        // Lyngso for the rest.
        for (int l = 6; l <= max_inter; ++l)
          mins[DP_P] = std::min(mins[DP_P], lyngso[st + 2][en - 2][l - 4] -
              g_internal_init[l - 4] + g_internal_init[l] + base_internal_loop);

        // Hairpin loops.
        mins[DP_P] = std::min(mins[DP_P], FastHairpin(st, en, hairpin_precomp));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = g_augubranch[stb][enb] + g_multiloop_hack_a;

        // (<   ><   >)
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + arr[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + arr[st + 2][en - 1][DP_U2] + g_dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        mins[DP_P] = std::min(mins[DP_P], base_branch_cost + arr[st + 1][en - 2][DP_U2] + g_dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        mins[DP_P] = std::min(mins[DP_P],
            base_branch_cost + arr[st + 2][en - 2][DP_U2] + g_terminal[stb][st1b][en1b][enb]);

        // (.(   ).   ) Left right coax
        for (auto cand : cand_st[CAND_P_MISMATCH])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy + arr[cand.idx + 1][en - 1][DP_U]);
        // (.(   )   .) Left outer coax
        auto outer_coax = energy::MismatchCoaxial(stb, st1b, en1b, enb);
        for (auto cand : cand_st[CAND_P_OUTER])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy - g_min_mismatch_coax +
              outer_coax + arr[cand.idx + 1][en - 2][DP_U]);
        // ((   )   ) Left flush coax
        for (auto cand : cand_st[CAND_P_FLUSH])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy - g_min_flush_coax +
              g_stack[stb][st1b][r[cand.idx]][enb] + arr[cand.idx + 1][en - 1][DP_U]);

        // (   .(   ).) Right left coax
        for (auto cand : p_cand_en[CAND_EN_P_MISMATCH][en])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy + arr[st + 1][cand.idx - 1][DP_U]);
        // (.   (   ).) Right outer coax
        for (auto cand : p_cand_en[CAND_EN_P_OUTER][en])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy - g_min_mismatch_coax +
              outer_coax + arr[st + 2][cand.idx - 1][DP_U]);
        // (   (   )) Right flush coax
        for (auto cand : p_cand_en[CAND_EN_P_FLUSH][en])
          mins[DP_P] = std::min(mins[DP_P], base_branch_cost + cand.energy - g_min_flush_coax +
              g_stack[stb][r[cand.idx]][en1b][enb] + arr[st + 1][cand.idx - 1][DP_U]);

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
        auto val = cand.energy - g_min_mismatch_coax + arr[cand.idx + 1][en][DP_U_RCOAX];
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
      // TODO refactor these comments.

      // Now build the candidates arrays based off the current pair [st, en].
      // In general, the idea is to see if there exists something we could replace a structure with that is as good.
      // e.g. we could replace (   )3' in unpaired with the equivalent U since we know the energy for (   )3' -
      // it is self contained. In some cases we use the minimum possible energy if we don't know the energy exactly
      // for a structure (e.g. RCOAX).
      // These orderings are useful to remember:
      // U <= U_WC, U_GU, U2

      energy_t cand_st_mins[] = {MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E, MAX_E};
      static_assert(sizeof(cand_st_mins) / sizeof(cand_st_mins[0]) == CAND_SIZE, "array wrong size");

      // Unpaired cases. These store the best pairs u_cand that begin at st.
      // begin means that the whole interaction starts at st. e.g. .(   ). starts one before the paren.
      // (   ) - Normal - U, U2
      auto normal_base = arr[st][en][DP_P] + g_augubranch[stb][enb];
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

      // TODO put if statments around these to make them only execute when actually possible?
      // TODO do assignemnt of arr[st][en] inside these if statements - or use current best; don't push same thing on
      // TODO order these if statements so hte most liekly to be strong go first.
      // multiple times? when this is split into multiple candidate lists still can have minimum over all of them for
      // monotonicity``
      // Can only merge candidate lists for monotonicity if the right part of the pivot is the same (from the same array).
      // Can only apply monotonicity optimisation to ones ending with min(U, 0).
      // (   ). - 3' - U, U2
      auto dangle3_base = arr[st][en - 1][DP_P] + g_augubranch[stb][en1b] + g_dangle3[en1b][enb][stb];
      if (dangle3_base < arr[st][en][DP_U] && dangle3_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = dangle3_base;
      // .(   ) - 5' - U, U2
      auto dangle5_base = arr[st + 1][en][DP_P] + g_augubranch[st1b][enb] + g_dangle5[enb][stb][st1b];
      if (dangle5_base < arr[st][en][DP_U] && dangle5_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = dangle5_base;
      // .(   ). - Terminal mismatch - U, U2
      auto terminal_base = arr[st + 1][en - 1][DP_P] + g_augubranch[st1b][en1b] + g_terminal[en1b][enb][stb][st1b];
      if (terminal_base < arr[st][en][DP_U] && terminal_base < cand_st_mins[CAND_U])
        cand_st_mins[CAND_U] = terminal_base;
      // .(   ).<(   ) > - Left coax - U, U2
      auto lcoax_base = arr[st + 1][en - 1][DP_P] + g_augubranch[st1b][en1b] +
          energy::MismatchCoaxial(en1b, enb, stb, st1b);
      if (lcoax_base < arr[st][en][DP_U])
        cand_st[CAND_U_LCOAX].push_back({lcoax_base, en});
      // (   ).<(   ). > Right coax forward - U, U2
      // This is probably better than having four candidate lists for each possible mismatch (TODO check this).
      // TODO remember to subtract off g_coax_mismatch_non_contiguous when computing after saving value
      auto rcoaxf_base = arr[st][en - 1][DP_P] + g_augubranch[stb][en1b] + g_min_mismatch_coax;
      if (rcoaxf_base < arr[st][en][DP_U])
        cand_st[CAND_U_RCOAX_FWD].push_back({rcoaxf_base, en});

      // (   ).<( * ). > Right coax backward - RCOAX
      // Again, we can't replace RCOAX with U, we'd have to replace it with RCOAX, so compare to itself.
      if (st > 0) {
        auto rcoaxb_base = arr[st][en - 1][DP_P] + g_augubranch[stb][en1b] + energy::MismatchCoaxial(
            en1b, enb, r[st - 1], stb);
        if (rcoaxb_base < arr[st][en][DP_U_RCOAX] && rcoaxb_base < cand_st_mins[CAND_U_RCOAX])
          cand_st_mins[CAND_U_RCOAX] = rcoaxb_base;
        // Base case.
        arr[st][en][DP_U_RCOAX] = std::min(arr[st][en][DP_U_RCOAX], rcoaxb_base);
      }

      // (   )<(   ) > Flush coax - U, U2
      if (en + 1 < N) {
        auto enr1b = r[en + 1];
        auto wc_flush_base = arr[st][en][DP_P] + g_augubranch[stb][enb] + g_stack[enb][enr1b][enr1b ^ 3][stb];
        auto gu_flush_base = arr[st][en][DP_P] + g_augubranch[stb][enb] + g_stack[enb][enr1b][enr1b ^ 1][stb];
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
      auto plocoax_base = arr[st + 2][en][DP_P] + g_augubranch[st2b][enb] + g_min_mismatch_coax;
      if (plocoax_base < arr[st + 1][en][DP_U])
        cand_st[CAND_P_OUTER].push_back({plocoax_base, en});
      // (.   (   ).) Right outer coax
      auto procoax_base = arr[st][en - 2][DP_P] + g_augubranch[stb][en2b] + g_min_mismatch_coax;
      if (procoax_base < arr[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_OUTER][en].push_back({procoax_base, st});
      // (.(   ).   ) Left right coax
      auto plrcoax_base = arr[st + 2][en - 1][DP_P] + g_augubranch[st2b][en1b] +
          energy::MismatchCoaxial(en1b, enb, st1b, st2b);
      if (plrcoax_base < arr[st + 1][en][DP_U])
        cand_st[CAND_P_MISMATCH].push_back({plrcoax_base, en});
      // (   .(   ).) Right left coax
      auto prlcoax_base = arr[st + 1][en - 2][DP_P] + g_augubranch[st1b][en2b] +
          energy::MismatchCoaxial(en2b, en1b, stb, st1b);
      if (prlcoax_base < arr[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_MISMATCH][en].push_back({prlcoax_base, st});
      // ((   )   ) Left flush coax
      auto plfcoax_base = arr[st + 1][en][DP_P] + g_augubranch[st1b][enb] + g_min_flush_coax;
      if (plfcoax_base < arr[st + 1][en][DP_U])
        cand_st[CAND_P_FLUSH].push_back({plfcoax_base, en});
      // (   (   )) Right flush coax
      auto prfcoax_base = arr[st][en - 1][DP_P] + g_augubranch[stb][en1b] + g_min_flush_coax;
      if (prfcoax_base < arr[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_FLUSH][en].push_back({prfcoax_base, st});


      // Add potentials to the candidate lists.
      for (int i = 0; i < CAND_SIZE; ++i) {
        if (cand_st_mins[i] < constants::CAP_E &&
            (cand_st[i].empty() || cand_st_mins[i] < cand_st[i].back().energy))
          cand_st[i].push_back({cand_st_mins[i], en});
      }
    }
  }
  return arr;
}

}
}
