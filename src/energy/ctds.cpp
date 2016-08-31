#include <algorithm>
#include <memory>
#include <cmath>
#include <cstdio>
#include "energy/ctds.h"
#include "energy/energy.h"
#include "energy/energy_globals.h"
#include "constants.h"

namespace memerna {
namespace energy {

// Computes the optimal arrangement of coaxial stackings, terminal mismatches, and dangles (CTD).
// This DP needs to be run four times. The series of branches is actually cyclic, and there are two types of interactions
// that can occur between branches. Interactions between the branches themselves, and interactions between them and the
// bases next to them. N.B. interact is not used in the chemistry sense here.
// The DP needs to be run with the outer branch at the start and the end, and also with the unpaired base on the left
// of the first branch (if it exists) not used and potentially used. Running with the outer branch at both the start and
// the end allows coaxial stacking interactions between the outer branch and both adjacent branches. Running with having the
// left unpaired base both not used lets the last branch potentially consume it (wrapping around).
//
// The state holds which branches we're up to, |i|, and whether the unpaired base on the left was consumed (if it exists).
//
// Note that we can't form two dangles attached to the same stem; that's a terminal mismatch.

// unpaired base will never have been consumed by a 5' dangle.
#define UPDATE_CACHE(used, idx, cur_used, cur_idx, value, reason) \
  do { \
    energy_t macro_upd_value_ = (value); \
    if (cache[cur_used][cur_idx] + macro_upd_value_ < cache[used][idx]) { \
      cache[used][idx] = cache[cur_used][cur_idx] + macro_upd_value_; \
      back[used][idx] = std::make_tuple(cur_used, cur_idx, macro_upd_value_, reason); \
    } \
  } while (0)
// TODO remove outer_idx, do reversal in calling code
energy_t ComputeOptimalCtd(const std::deque<int>& branches,
    int outer_idx, bool use_first_lu, branch_ctd_t* ctd_energies) {
  int N = int(branches.size());
  int RSZ = int(r.size());
  assert(outer_idx == 0 || outer_idx == N - 1 || outer_idx == -1);
  // Could be on the exterior loop.
  assert(N >= 3 || outer_idx == -1);
  if (N < 1) return 0;

  // cache[used][i]
  std::vector<int> cache[2] = {
      std::vector<int>(size_t(N + 1), constants::MAX_E),
      std::vector<int>(size_t(N + 1), constants::MAX_E)
  };
  std::vector<std::tuple<bool, int, energy_t, Ctd>> back[2] = {
      std::vector<std::tuple<bool, int, energy_t, Ctd>>(size_t(N + 1), std::make_tuple(false, -1, 0, CTD_NA)),
      std::vector<std::tuple<bool, int, energy_t, Ctd>>(size_t(N + 1), std::make_tuple(false, -1, 0, CTD_NA))
  };

  cache[0][0] = cache[1][0] = 0;
  int first_lui = -1, last_rui = -1;
  if (outer_idx == 0) {
    first_lui = p[branches[0]] - 1;
    last_rui = p[branches[N - 1]] + 1;
  } else if (outer_idx == N - 1) {
    first_lui = branches[0] - 1;
    last_rui = branches[N - 1] + 1;
  }

  // Precompute data about the unpaired bases.
  std::vector<int> li((size_t(N))), ri((size_t(N))), lui((size_t(N))), rui((size_t(N)));
  std::vector<bool> lu_exists((size_t(N))), lu_usable((size_t(N))),
      ru_exists((size_t(N))), ru_usable((size_t(N))), ru_shared((size_t(N)));
  for (int i = 0; i < N; ++i) {
    li[i] = branches[i];
    ri[i] = p[branches[i]];
    assert(ri[i] != -1);
    if (i == outer_idx)
      std::swap(li[i], ri[i]);
    lui[i] = li[i] - 1;
    rui[i] = ri[i] + 1;
    // If |use_first_lu|, then if the left unpaired base is the same as the last branch's right unpaired base,
    // then we can't use it (as it could be used at the end by a terminal mismatch, dangle, right facing coaxial stack,
    // etc). This is because the loop is cyclic.
    lu_exists[i] = lui[i] >= 0 && lui[i] < RSZ && p[lui[i]] == -1;
    lu_usable[i] = lu_exists[i] && (lui[i] != last_rui || use_first_lu);
    ru_exists[i] = rui[i] >= 0 && rui[i] < RSZ && p[rui[i]] == -1;
    ru_usable[i] = ru_exists[i] && (rui[i] != first_lui || !use_first_lu);
    ru_shared[i] = ru_exists[i] && rui[i] < RSZ - 1 && p[rui[i] + 1] != -1;
  }

  for (int i = 0; i < N; ++i) {
    base_t lb = r[li[i]], rb = r[ri[i]], lub = -1, rub = -1;
    if (lu_exists[i]) lub = r[lui[i]];
    if (ru_exists[i]) rub = r[rui[i]];

    // Flush coaxial stacking. Requires that ru not exist (i.e. adjacent branches) and this not be the last branch.
    if (!ru_exists[i] && i != N - 1) {
      energy_t coax = g_stack[rb][r[li[i + 1]]][r[ri[i + 1]]][lb];
      // When the next branch is consumed by this coaxial stack, it can no longer interact with anything, so
      // just skip to i + 2.
      UPDATE_CACHE(0, i + 2, 0, i, coax, CTD_FLUSH_COAX_WITH_NEXT);
      if (lu_exists[i]) {
        // If lu exists, and it was used, then it's fine to coaxially stack. If |used| were true but lu didn't exist then
        // we couldn't coaxially stack as the current branch would already have been involved in one, though.
        UPDATE_CACHE(0, i + 2, 1, i, coax, CTD_FLUSH_COAX_WITH_NEXT);
      }
    }

    if (lu_usable[i] && ru_usable[i]) {
      // Terminal mismatch, requires lu_exists, ru_exists, and that we didn't use left.
      // Consumes ru, so if it was shared, use it.
      UPDATE_CACHE(ru_shared[i], i + 1, 0, i, g_terminal[rb][rub][lub][lb], CTD_TERMINAL_MISMATCH);

      // Mismatch mediated coaxial stacking, left facing (uses the branch we're currently at).
      // Requires lu_usable, ru_usable, ru_shared, and left not used. Consumes ru.
      // Skip to the branch after next since the next branch can't be involved in any more interactions anyway:
      // its left pair is consumed, and its right pair can't dangle towards it.
      if (ru_shared[i] && i != N - 1) {
        energy_t left_coax = energy::MismatchCoaxial(rb, rub, lub, lb);
        UPDATE_CACHE(0, i + 2, 0, i, left_coax, CTD_LEFT_MISMATCH_COAX_WITH_NEXT);
      }
    }

    // Right consuming cases.
    if (ru_usable[i]) {
      // Right dangle (3').
      // Only requires ru_exists so handle where left is both used and not used.
      UPDATE_CACHE(ru_shared[i], i + 1, 0, i, g_dangle3[rb][rub][lb], CTD_3_DANGLE);
      UPDATE_CACHE(ru_shared[i], i + 1, 1, i, g_dangle3[rb][rub][lb], CTD_3_DANGLE);

      // Mismatch mediated coaxial stacking, right facing (uses the next branch).
      // Requires ru_exists, ru_shared. Consumes ru and rru.
      if (ru_shared[i] && i != N - 1 && ru_usable[i + 1]) {
        energy_t right_coax = energy::MismatchCoaxial(r[ri[i + 1]], r[rui[i + 1]], rub, r[li[i + 1]]);

        UPDATE_CACHE(ru_shared[i + 1], i + 2, 0, i, right_coax, CTD_RIGHT_MISMATCH_COAX_WITH_NEXT);
        if (lu_exists[i]) {
          // In the case that lu doesn't exist but it is "used" it means this branch was consumed by a coaxial interaction
          // so don't use it.
          UPDATE_CACHE(ru_shared[i + 1], i + 2, 1, i, right_coax, CTD_RIGHT_MISMATCH_COAX_WITH_NEXT);
        }
      }
    }

    if (lu_usable[i]) {
      // 5' dangle.
      UPDATE_CACHE(0, i + 1, 0, i, g_dangle5[rb][lub][lb], CTD_5_DANGLE);
    }

    // Have the option of doing nothing.
    UPDATE_CACHE(0, i + 1, 0, i, 0, CTD_UNUSED);
    UPDATE_CACHE(0, i + 1, 1, i, 0, CTD_UNUSED);
  }

  std::tuple<bool, int, energy_t, Ctd> state{false, N, 0, CTD_NA};
  if (cache[1][N] < cache[0][N])
    state = std::make_tuple(true, N, 0, CTD_NA);
  // First state contains no real info, so go ahead one.
  state = back[std::get<0>(state)][std::get<1>(state)];
  assert(ctd_energies->empty());
  while (1) {
    bool used;
    int idx;
    energy_t energy;
    Ctd reason;
    std::tie(used, idx, energy, reason) = std::move(state);
    if (idx == -1) break;
    // We can skip a branch on coaxial stacking interactions, so make sure to insert the ctd energy for the branch.
    // We build ctd_energies backwards, so we need to insert this later branch first.
    // To get the PREV versions, just add one.
    if (reason == CTD_LEFT_MISMATCH_COAX_WITH_NEXT ||
        reason == CTD_RIGHT_MISMATCH_COAX_WITH_NEXT ||
        reason == CTD_FLUSH_COAX_WITH_NEXT)
      ctd_energies->emplace_back(Ctd(reason + 1), energy);
    ctd_energies->emplace_back(reason, energy);
    state = back[used][idx];
  }
  std::reverse(ctd_energies->begin(), ctd_energies->end());

  return std::min(cache[0][N], cache[1][N]);
}

#undef UPDATE_CACHE

void WriteCtdsForBranches(const std::deque<int>& branches, const branch_ctd_t& ctd_energies) {
  assert(branches.size() == ctd_energies.size());
  for (int i = 0; i < int(branches.size()); ++i) {
    ctds[branches[i]] = ctd_energies[i].first;
    ctds[p[branches[i]]] = ctd_energies[i].first;
  }
}

energy_t ReadCtdsForBranches(const std::deque<int>& branches, branch_ctd_t* ctd_energies) {
  energy_t total_energy = 0;
  int N = int(r.size());
  for (int i = 0; i < int(branches.size()); ++i) {
    int branch = branches[i];
    int prev_branch = i > 0 ? branches[i - 1] : -1;
    assert(ctds[branch] == ctds[p[branch]]);
    energy_t energy = 0;
    auto stb = r[branch], enb = r[p[branch]];
    // TODO swap if outer branch
    switch (ctds[branch]) {
      case CTD_UNUSED:
        break;
      case CTD_3_DANGLE:
        assert(branch > 0);
        energy = g_dangle3[enb][r[branch - 1]][stb];
        break;
      case CTD_5_DANGLE:
        assert(branch > 0);
        energy = g_dangle5[enb][r[branch - 1]][stb];
        break;
      case CTD_TERMINAL_MISMATCH:
        assert(p[branch] + 1 < N && branch > 0);
        energy = g_terminal[enb][r[p[branch] + 1]][r[branch - 1]][stb];
        break;
      case CTD_LEFT_MISMATCH_COAX_WITH_PREV:
        assert(prev_branch != -1);
        // .(   ).(   )
        energy = MismatchCoaxial(r[p[prev_branch]], r[p[prev_branch] + 1], r[prev_branch - 1], r[prev_branch]);
        ctd_energies->emplace_back(CTD_LEFT_MISMATCH_COAX_WITH_NEXT, energy);
        break;
      case CTD_RIGHT_MISMATCH_COAX_WITH_PREV:
        assert(prev_branch != -1 && branch > 0 && p[branch] + 1 < N);
        // (   ).(   ). or (.(   ).   )
        energy = MismatchCoaxial(enb, r[p[branch] + 1], r[branch - 1], stb);
        ctd_energies->emplace_back(CTD_RIGHT_MISMATCH_COAX_WITH_NEXT, energy);
        break;
      case CTD_FLUSH_COAX_WITH_PREV:
        assert(prev_branch != -1);
        // (   )(   ) or ((   )   ) or (   (   ))
        energy = g_stack[r[p[prev_branch]]][stb][enb][r[prev_branch]];
        break;
        // All these cases will be handled in the next branch (PREV).
      case CTD_LEFT_MISMATCH_COAX_WITH_NEXT:
      case CTD_RIGHT_MISMATCH_COAX_WITH_NEXT:
      case CTD_FLUSH_COAX_WITH_NEXT:
        break;
      default:
        verify_expr(false, "bug");  // Should never happen
    }
    ctd_energies->emplace_back(ctds[branch], energy);
  }
  return total_energy;
}

}
}
