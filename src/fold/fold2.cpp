#include "fold2.h"

namespace memerna {
namespace fold {

// Non continuous -2.1, -4 for WC, -16 for terminal mismatch.
const energy_t MIN_MISMATCH_COAX = -21 - 4 - 16;
const energy_t MIN_FLUSH_COAX = -34;

// Split candidates up into several lists.
// In general, for each array we need a new candidate list (except for U and U2 which mirror each other very
// closely). We also need another candidate list for forward RCOAX since we can't use its energy value directly, not
// knowing it. Same with flush coaxial stacks.
enum {
  CAND_P_MISMATCH,
  CAND_P_OUTER,
  CAND_P_FLUSH,
  CAND_U,
  CAND_U_LCOAX,
  CAND_U_RCOAX_FWD,
  CAND_U_FLUSH,
  CAND_U_WC,
  CAND_U_GU,
  CAND_U_RCOAX,
  CAND_SIZE
};

enum {
  CAND_EN_P_MISMATCH,
  CAND_EN_P_OUTER,
  CAND_EN_P_FLUSH,
  CAND_EN_SIZE
};

struct cand_t {
  energy_t energy;
  int idx;
};

array3d_t<energy_t, DP_SIZE> ComputeTables2() {
  int N = int(r.size());
  // Automatically initialised to MAX_E.
  array3d_t<energy_t, DP_SIZE> arr(r.size() + 1);
  // For now split up cases into paired and unpaired candidates. For more optimisations, split this up more later.
  // In the full optimisation, these are strictly decreasing (since 0 unpaired base cost).
  // In the full optimisation, also include the computed values (except for RCOAX which is not possible).
  std::vector<std::vector<cand_t>> p_cand_en[CAND_EN_SIZE];
  for (auto& i : p_cand_en) i.resize(r.size());
  std::vector<cand_t> cand_st[CAND_SIZE];

  static_assert(constants::HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int st = N - 1; st >= 0; --st) {
    for (int i = 0; i < CAND_SIZE; ++i)
      cand_st[i].clear();
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      base_t stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];

      // Update paired - only if can actually pair.
      energy_t p_min = constants::MAX_E;
      if (CanPair(r[st], r[en]) && IsNotLonely(st, en)) {
        int max_inter = std::min(constants::TWOLOOP_MAX_SZ, en - st - constants::HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (arr[ist][ien][DP_P] < constants::CAP_E)
              p_min = std::min(p_min, FastTwoLoop(st, en, ist, ien) + arr[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        p_min = std::min(p_min, energy::HairpinEnergy(st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = AUGUBRANCH[stb][enb] + MULTILOOP_A;

        // (<   ><   >)
        p_min = std::min(p_min, base_branch_cost + arr[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        p_min = std::min(p_min, base_branch_cost + arr[st + 2][en - 1][DP_U2] + dangle3_e[stb][st1b][enb]);
        // (<   ><   >5) 5'
        p_min = std::min(p_min, base_branch_cost + arr[st + 1][en - 2][DP_U2] + dangle5_e[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        p_min = std::min(p_min, base_branch_cost + arr[st + 2][en - 2][DP_U2] + terminal_e[stb][st1b][en1b][enb]);

        // (.(   ).   ) Left right coax
        for (auto cand : cand_st[CAND_P_MISMATCH])
          p_min = std::min(p_min, base_branch_cost + cand.energy + arr[cand.idx + 1][en - 1][DP_U]);
        // (.(   )   .) Left outer coax
        auto outer_coax = energy::MismatchMediatedCoaxialEnergy(stb, st1b, en1b, enb);
        for (auto cand : cand_st[CAND_P_OUTER])
          p_min = std::min(p_min, base_branch_cost + cand.energy - MIN_MISMATCH_COAX +
              outer_coax + arr[cand.idx + 1][en - 2][DP_U]);
        // ((   )   ) Left flush coax
        for (auto cand : cand_st[CAND_P_FLUSH])
          p_min = std::min(p_min, base_branch_cost + cand.energy - MIN_FLUSH_COAX +
              stacking_e[stb][st1b][r[cand.idx]][enb] + arr[cand.idx + 1][en - 1][DP_U]);

        // (   .(   ).) Right left coax
        for (auto cand : p_cand_en[CAND_EN_P_MISMATCH][en])
          p_min = std::min(p_min, base_branch_cost + cand.energy + arr[st + 1][cand.idx - 1][DP_U]);
        // (.   (   ).) Right outer coax
        for (auto cand : p_cand_en[CAND_EN_P_OUTER][en])
          p_min = std::min(p_min, base_branch_cost + cand.energy - MIN_MISMATCH_COAX +
              outer_coax + arr[st + 2][cand.idx - 1][DP_U]);
        // (   (   )) Right flush coax
        for (auto cand : p_cand_en[CAND_EN_P_FLUSH][en])
          p_min = std::min(p_min, base_branch_cost + cand.energy - MIN_FLUSH_COAX +
              stacking_e[stb][r[cand.idx]][en1b][enb] + arr[st + 1][cand.idx - 1][DP_U]);

        arr[st][en][DP_P] = p_min;
      }
      energy_t u_min = constants::MAX_E, u2_min = constants::MAX_E,
          u_gu_min = constants::MAX_E, u_wc_min = constants::MAX_E, u_rcoax_min = constants::MAX_E;
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u_min = std::min(u_min, arr[st + 1][en][DP_U]);
        u2_min = std::min(u2_min, arr[st + 1][en][DP_U2]);
      }
      for (auto cand : cand_st[CAND_U]) {
        u_min = std::min(u_min, cand.energy + std::min(arr[cand.idx + 1][en][DP_U], 0));
        u2_min = std::min(u2_min, cand.energy + arr[cand.idx + 1][en][DP_U]);
      }
      for (auto cand : cand_st[CAND_U_LCOAX]) {
        auto val = cand.energy + std::min(arr[cand.idx + 1][en][DP_U_WC], arr[cand.idx + 1][en][DP_U_GU]);
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
      }
      for (auto cand : cand_st[CAND_U_RCOAX_FWD]) {
        auto val = cand.energy - MIN_MISMATCH_COAX + arr[cand.idx + 1][en][DP_U_RCOAX];
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
      }
      for (auto cand : cand_st[CAND_U_FLUSH]) {
        // (   )<(   ) > Flush coax - U
        // stb piv piv + 1
        if (cand.idx + 1 < en) {
          auto pb = r[cand.idx], pr1b = r[cand.idx + 1];
          auto val =
              cand.energy - MIN_FLUSH_COAX + stacking_e[pb][pr1b][pr1b ^ 3][stb] + arr[cand.idx + 1][en][DP_U_WC];
          u_min = std::min(u_min, val);
          u2_min = std::min(u2_min, val);
          if (pr1b == G || pr1b == U) {
            val = cand.energy - MIN_FLUSH_COAX + stacking_e[pb][pr1b][pr1b ^ 1][stb] + arr[cand.idx + 1][en][DP_U_GU];
            u_min = std::min(u_min, val);
            u2_min = std::min(u2_min, val);
          }
        }
      }
      for (auto cand : cand_st[CAND_U_WC])
        u_wc_min = std::min(u_wc_min, cand.energy + std::min(arr[cand.idx + 1][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_GU])
        u_gu_min = std::min(u_gu_min, cand.energy + std::min(arr[cand.idx + 1][en][DP_U], 0));
      for (auto cand : cand_st[CAND_U_RCOAX]) {
        // (   ).<( * ). > Right coax backward
        assert(st > 0);
        u_rcoax_min = std::min(u_rcoax_min, cand.energy + std::min(arr[cand.idx + 1][en][DP_U], 0));
      }

      // Set these so we can use sparse folding.
      arr[st][en][DP_U] = u_min;
      arr[st][en][DP_U2] = u2_min;
      arr[st][en][DP_U_WC] = u_wc_min;
      arr[st][en][DP_U_GU] = u_gu_min;
      arr[st][en][DP_U_RCOAX] = u_rcoax_min;
      // TODO refactor these comments.

      // Now build the candidates arrays based off the current pair [st, en].
      // When we build the unpaired and paired arrays, we assume there is a pair starting at st, and then
      // try every pivot. Since U[st][en] has to be <= than P[st][en] (since U it has the option of pairing there
      // anyway), we only want to try using P[st][en] when it's better than U[st][en]. Otherwise we will have already
      // tried a P[st][en'] with en' < en that U[st][en] covers, which is better.
      // We need to consider each case covering the whole of st, en as part of the base cases for this to work.
      // This is only for interactions that can exist by themselves (i.e. not coaxial stacking).
      // TODO check that this correctly does no work for impossible pairings!!

      // In general, the idea is to see if there exists something we could replace a structure with that is as good.
      // e.g. we could replace (   )3' in unpaired with the equivalent U since we know the energy for (   )3' -
      // it is self contained. In some cases we use the minimum possible energy if we don't know the energy exactly
      // for a structure (e.g. RCOAX).

      // "Replaceability"

      // These orderings are useful to remember:
      // U <= U_WC, U_GU, U2

      // Unpaired cases. These store the best pairs u_cand that begin at st.
      // begin means that the whole interaction starts at st. e.g. .(   ). starts one before the paren.
      // This means the current pair [st, en] can update locations at st - 1.
      // (   ) - Normal - U, U2, U_WC?, U_GU?
      auto normal_base = arr[st][en][DP_P] + AUGUBRANCH[stb][enb];
      if (normal_base < arr[st][en][DP_U])
        cand_st[CAND_U].push_back({normal_base, en});

      // For U_GU and U_WC, they can't be replaced with DP_U, so we need to compare them to something they can be
      // replaced with, i.e. themselves.
      if (IsGu(stb, enb)) {
        if (normal_base < arr[st][en][DP_U_GU])
          cand_st[CAND_U_GU].push_back({normal_base, en});
        // Base case.
        arr[st][en][DP_U_GU] = std::min(arr[st][en][DP_U_GU], normal_base);
      } else {
        if (normal_base < arr[st][en][DP_U_WC])
          cand_st[CAND_U_WC].push_back({normal_base, en});
        // Base case.
        arr[st][en][DP_U_WC] = std::min(arr[st][en][DP_U_WC], normal_base);
      }

      // TODO put if statments around these to make them only execute when actually possible.
      // TODO save these values so don't have to recalc (except for RCOAX stack - can just do subtraction and add)
      // TODO do assignemnt of arr[st][en] inside these if statements - or use current best; don't push same thing on
      // multiple times? when this is split into multiple candidate lists still can have minimum over all of them for
      // monotonicity``
      // Can only merge candidate lists for monotonicity if the right part of the pivot is the same (from the same array).
      // (   ). - 3' - U, U2
      auto dangle3_base = arr[st][en - 1][DP_P] + AUGUBRANCH[stb][en1b] + dangle3_e[en1b][enb][stb];
      if (dangle3_base < arr[st][en][DP_U])
        cand_st[CAND_U].push_back({dangle3_base, en});
      // .(   ) - 5' - U, U2
      auto dangle5_base = arr[st + 1][en][DP_P] + AUGUBRANCH[st1b][enb] + dangle5_e[enb][stb][st1b];
      if (dangle5_base < arr[st][en][DP_U])
        cand_st[CAND_U].push_back({dangle5_base, en});
      // .(   ). - Terminal mismatch - U, U2
      auto terminal_base = arr[st + 1][en - 1][DP_P] + AUGUBRANCH[st1b][en1b] + terminal_e[en1b][enb][stb][st1b];
      if (terminal_base < arr[st][en][DP_U])
        cand_st[CAND_U].push_back({terminal_base, en});
      // .(   ).<(   ) > - Left coax - U, U2
      auto lcoax_base = arr[st + 1][en - 1][DP_P] + AUGUBRANCH[st1b][en1b] +
          energy::MismatchMediatedCoaxialEnergy(en1b, enb, stb, st1b);
      if (lcoax_base < arr[st][en][DP_U])
        cand_st[CAND_U_LCOAX].push_back({lcoax_base, en});
      // (   ).<(   ). > Right coax forward - U, U2
      // This is probably better than having four candidate lists for each possible mismatch (TODO check this).
      // TODO remember to subtract off coax_mismatch_non_contiguous when computing after saving value
      auto rcoaxf_base = arr[st][en - 1][DP_P] + AUGUBRANCH[stb][en1b] + MIN_MISMATCH_COAX;
      if (rcoaxf_base < arr[st][en][DP_U])
        cand_st[CAND_U_RCOAX_FWD].push_back({rcoaxf_base, en});

      // (   ).<( * ). > Right coax backward - RCOAX
      // Again, we can't replace RCOAX with U, we'd have to replace it with RCOAX, so compare to itself.
      if (st > 0) {
        auto rcoaxb_base = arr[st][en - 1][DP_P] + AUGUBRANCH[stb][en1b] + energy::MismatchMediatedCoaxialEnergy(
            en1b, enb, r[st - 1], stb);
        if (rcoaxb_base < arr[st][en][DP_U_RCOAX])
          cand_st[CAND_U_RCOAX].push_back({rcoaxb_base, en});
        // Base case.
        arr[st][en][DP_U_RCOAX] = std::min(arr[st][en][DP_U_RCOAX], rcoaxb_base);
      }

      // (   )<(   ) > Flush coax - U, U2
      // We could take the min of if it were a Watson-Crick or GU pair for stacking, but then we would have to
      // be very careful when keeping this candidate list monotonic, since stacks could have less or more energy
      // than we expect.
      auto flush_base = arr[st][en][DP_P] + AUGUBRANCH[stb][enb] + MIN_FLUSH_COAX;
      if (flush_base < arr[st][en][DP_U])
        cand_st[CAND_U_FLUSH].push_back({flush_base, en});

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
      auto plocoax_base = arr[st + 2][en][DP_P] + AUGUBRANCH[st2b][enb] + MIN_MISMATCH_COAX;
      if (plocoax_base < arr[st + 1][en][DP_U])
        cand_st[CAND_P_OUTER].push_back({plocoax_base, en});
      // (.   (   ).) Right outer coax
      auto procoax_base = arr[st][en - 2][DP_P] + AUGUBRANCH[stb][en2b] + MIN_MISMATCH_COAX;
      if (procoax_base < arr[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_OUTER][en].push_back({procoax_base, st});
      // (.(   ).   ) Left right coax
      auto plrcoax_base = arr[st + 2][en - 1][DP_P] + AUGUBRANCH[st2b][en1b] +
          energy::MismatchMediatedCoaxialEnergy(en1b, enb, st1b, st2b);
      if (plrcoax_base < arr[st + 1][en][DP_U])
        cand_st[CAND_P_MISMATCH].push_back({plrcoax_base, en});
      // (   .(   ).) Right left coax
      auto prlcoax_base = arr[st + 1][en - 2][DP_P] + AUGUBRANCH[st1b][en2b] +
          energy::MismatchMediatedCoaxialEnergy(en2b, en1b, stb, st1b);
      if (prlcoax_base < arr[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_MISMATCH][en].push_back({prlcoax_base, st});
      // ((   )   ) Left flush coax
      auto plfcoax_base = arr[st + 1][en][DP_P] + AUGUBRANCH[st1b][enb] + MIN_FLUSH_COAX;
      if (plfcoax_base < arr[st + 1][en][DP_U])
        cand_st[CAND_P_FLUSH].push_back({plfcoax_base, en});
      // (   (   )) Right flush coax
      auto prfcoax_base = arr[st][en - 1][DP_P] + AUGUBRANCH[stb][en1b] + MIN_FLUSH_COAX;
      if (prfcoax_base < arr[st][en - 1][DP_U])
        p_cand_en[CAND_EN_P_FLUSH][en].push_back({prfcoax_base, st});
    }
  }
  return arr;
}

}
}
