// Copyright 2016 Eliot Courtney.
#include "energy/energy_internal.h"

#include <algorithm>

#include "parsing.h"

namespace memerna {
namespace energy {
namespace internal {

// Computes the optimal arrangement of coaxial stackings, terminal mismatches, and dangles (CTD).
// This DP needs to be run four times. The series of branches is actually cyclic, and there are two
// types of interactions that can occur between branches. Interactions between the branches
// themselves, and interactions between them and the bases next to them. N.B. interact is not
// used in the chemistry sense here. The DP needs to be run with the outer branch at the start
// and the end, and also with the unpaired base on the left of the first branch (if it exists)
// not used and potentially used. Running with the outer branchat both the start and
// the end allows coaxial stacking interactions between the outer branch and both adjacent branches.
// Running with having the left unpaired base both not used lets the last branch potentially
// consume it (wrapping around).
//
// The state holds which branches we're up to, |i|, and whether the unpaired base on the left was
// consumed (if it exists).
//
// Note that we can't form two dangles attached to the same stem; that's a terminal mismatch.

// unpaired base will never have been consumed by a 5' dangle.
#define UPDATE_CACHE(used, idx, cur_used, cur_idx, value, reason)                     \
  do {                                                                                \
    energy_t macro_upd_value_ = (value);                                              \
    if (cache[cur_used][cur_idx] + macro_upd_value_ < cache[used][idx]) {             \
      cache[used][idx] = cache[cur_used][cur_idx] + macro_upd_value_;                 \
      back[used][idx] = std::make_tuple(cur_used, cur_idx, macro_upd_value_, reason); \
    }                                                                                 \
  } while (0)
energy_t ComputeOptimalCtd(const secondary_t& secondary, const EnergyModel& em,
    const std::deque<int>& branches, bool use_first_lu, branch_ctd_t& branch_ctds) {
  const auto& r = secondary.r;
  const auto& p = secondary.p;
  int N = int(branches.size());
  int RSZ = int(r.size());
  assert(branch_ctds.empty());
  // Could be on the exterior loop with a branch (0, N - 1).
  if (N < 1) return 0;

  // cache[used][i]
  std::vector<int> cache[2] = {
      std::vector<int>(size_t(N + 1), MAX_E), std::vector<int>(size_t(N + 1), MAX_E)};
  std::vector<std::tuple<bool, int, energy_t, Ctd>> back[2] = {
      std::vector<std::tuple<bool, int, energy_t, Ctd>>(
          size_t(N + 1), std::make_tuple(false, -1, 0, CTD_NA)),
      std::vector<std::tuple<bool, int, energy_t, Ctd>>(
          size_t(N + 1), std::make_tuple(false, -1, 0, CTD_NA))};

  cache[0][0] = cache[1][0] = 0;
  int first_lui = branches[0] - 1, last_rui = p[branches[N - 1]] + 1;

  // Precompute data about the unpaired bases.
  std::vector<int> li((size_t(N))), ri((size_t(N))), lui((size_t(N))), rui((size_t(N)));
  std::vector<bool> lu_exists((size_t(N))), lu_usable((size_t(N))), ru_exists((size_t(N))),
      ru_usable((size_t(N))), ru_shared((size_t(N)));
  for (int i = 0; i < N; ++i) {
    li[i] = branches[i];
    ri[i] = p[branches[i]];
    assert(ri[i] != -1);
    lui[i] = li[i] - 1;
    rui[i] = ri[i] + 1;
    // If |use_first_lu|, then if the left unpaired base is the same as the last branch's right
    // unpaired base, then we can't use it (as it could be used at the end by a terminal mismatch,
    // dangle, right facing coaxial stack, etc). This is because the loop is cyclic.
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

    // Flush coaxial stacking. Requires that ru not exist (i.e. adjacent branches) and this not be
    // the last branch.
    if (!ru_exists[i] && i != N - 1) {
      energy_t coax = em.stack[rb][r[li[i + 1]]][r[ri[i + 1]]][lb];
      // When the next branch is consumed by this coaxial stack, it can no longer interact with
      // anything, so just skip to i + 2.
      UPDATE_CACHE(0, i + 2, 0, i, coax, CTD_FCOAX_WITH_NEXT);
      if (lu_exists[i]) {
        // If lu exists, and it was used, then it's fine to coaxially stack. If |used| were true but
        // lu didn't exist then we couldn't coaxially stack as the current branch would already
        // have been involved in one, though.
        UPDATE_CACHE(0, i + 2, 1, i, coax, CTD_FCOAX_WITH_NEXT);
      }
    }

    if (lu_usable[i] && ru_usable[i]) {
      // Terminal mismatch, requires lu_exists, ru_exists, and that we didn't use left.
      // Consumes ru, so if it was shared, use it.
      UPDATE_CACHE(ru_shared[i], i + 1, 0, i, em.terminal[rb][rub][lub][lb], CTD_MISMATCH);

      // Mismatch mediated coaxial stacking, left facing (uses the branch we're currently at).
      // Requires lu_usable, ru_usable, ru_shared, and left not used. Consumes ru.
      // Skip to the branch after next since the next branch can't be involved in any more
      // interactions anyway: its left pair is consumed, and its right pair can't dangle towards it.
      if (ru_shared[i] && i != N - 1) {
        energy_t left_coax = em.MismatchCoaxial(rb, rub, lub, lb);
        UPDATE_CACHE(0, i + 2, 0, i, left_coax, CTD_LCOAX_WITH_NEXT);
      }
    }

    // Right consuming cases.
    if (ru_usable[i]) {
      // Right dangle (3').
      // Only requires ru_exists so handle where left is both used and not used.
      UPDATE_CACHE(ru_shared[i], i + 1, 0, i, em.dangle3[rb][rub][lb], CTD_3_DANGLE);
      UPDATE_CACHE(ru_shared[i], i + 1, 1, i, em.dangle3[rb][rub][lb], CTD_3_DANGLE);

      // Mismatch mediated coaxial stacking, right facing (uses the next branch).
      // Requires ru_exists, ru_shared. Consumes ru and rru.
      if (ru_shared[i] && i != N - 1 && ru_usable[i + 1]) {
        energy_t right_coax = em.MismatchCoaxial(r[ri[i + 1]], r[rui[i + 1]], rub, r[li[i + 1]]);

        UPDATE_CACHE(ru_shared[i + 1], i + 2, 0, i, right_coax, CTD_RCOAX_WITH_NEXT);
        if (lu_exists[i]) {
          // In the case that lu doesn't exist but it is "used" it means this branch was consumed by
          // a coaxial interaction so don't use it.
          UPDATE_CACHE(ru_shared[i + 1], i + 2, 1, i, right_coax, CTD_RCOAX_WITH_NEXT);
        }
      }
    }

    if (lu_usable[i]) {
      // 5' dangle.
      UPDATE_CACHE(0, i + 1, 0, i, em.dangle5[rb][lub][lb], CTD_5_DANGLE);
    }

    // Have the option of doing nothing.
    UPDATE_CACHE(0, i + 1, 0, i, 0, CTD_UNUSED);
    UPDATE_CACHE(0, i + 1, 1, i, 0, CTD_UNUSED);
  }

  std::tuple<bool, int, energy_t, Ctd> state{false, N, 0, CTD_NA};
  if (cache[1][N] < cache[0][N]) state = std::make_tuple(true, N, 0, CTD_NA);
  // First state contains no real info, so go ahead one.
  state = back[std::get<0>(state)][std::get<1>(state)];
  assert(branch_ctds.empty());
  while (1) {
    bool used;
    int idx;
    energy_t energy;
    Ctd reason;
    std::tie(used, idx, energy, reason) = std::move(state);
    if (idx == -1) break;
    // We can skip a branch on coaxial stacking interactions, so make sure to insert the ctd energy
    // for the branch. We build branch_ctds backwards, so we need to insert this later branch first.
    // To get the PREV versions, just add one.
    if (reason == CTD_LCOAX_WITH_NEXT || reason == CTD_RCOAX_WITH_NEXT ||
        reason == CTD_FCOAX_WITH_NEXT)
      branch_ctds.emplace_back(Ctd(reason + 1), energy);
    branch_ctds.emplace_back(reason, energy);
    state = back[used][idx];
  }
  std::reverse(branch_ctds.begin(), branch_ctds.end());

  return std::min(cache[0][N], cache[1][N]);
}

#undef UPDATE_CACHE

void AddBranchCtdsToComputed(
    computed_t& computed, const std::deque<int>& branches, const branch_ctd_t& branch_ctds) {
  assert(branches.size() == branch_ctds.size());
  for (int i = 0; i < int(branches.size()); ++i) {
    assert(int(computed.base_ctds.size()) > branches[i] &&
        int(computed.base_ctds.size()) > computed.s.p[branches[i]]);
    // Only write it into one side. If it's for an outer loop, it will be the right side, since we
    // swap the indices in that case.
    computed.base_ctds[branches[i]] = branch_ctds[i].first;
  }
}

energy_t GetBranchCtdsFromComputed(const computed_t& computed, const EnergyModel& em,
    const std::deque<int>& branches, branch_ctd_t& branch_ctds) {
  const auto& r = computed.s.r;
  const auto& p = computed.s.p;
  assert(branch_ctds.empty());
  energy_t total_energy = 0;
  // If we have an outer loop in |branches|, it is possible the first could refer to PREV, or the
  // last, to NEXT. In this case, we need to fix the branch_ctds so that the corresponding branch
  // ctd is on the right side. e.g. if the first element refers to PREV, we would put something
  // before it, but that actually needs to be at the end.
  bool rot_left = false;
  for (int i = 0; i < int(branches.size()); ++i) {
    const int branch = branches[i];
    const int prev_branch = i > 0 ? branches[i - 1] : branches.back();
    energy_t energy = 0;
    const auto stb = r[branch], enb = r[p[branch]];
    switch (computed.base_ctds[branch]) {
    case CTD_UNUSED: break;
    case CTD_3_DANGLE:
      assert(p[branch] + 1 < int(r.size()));
      energy = em.dangle3[enb][r[p[branch] + 1]][stb];
      break;
    case CTD_5_DANGLE:
      assert(branch > 0);
      energy = em.dangle5[enb][r[branch - 1]][stb];
      break;
    case CTD_MISMATCH:
      assert(p[branch] + 1 < int(r.size()) && branch > 0);
      energy = em.terminal[enb][r[p[branch] + 1]][r[branch - 1]][stb];
      break;
    case CTD_LCOAX_WITH_PREV:
      // .(   ).(   )
      assert(p[prev_branch] + 1 < int(r.size()) && prev_branch - 1 >= 0);
      energy = em.MismatchCoaxial(
          r[p[prev_branch]], r[p[prev_branch] + 1], r[prev_branch - 1], r[prev_branch]);
      branch_ctds.emplace_back(CTD_LCOAX_WITH_NEXT, energy);
      rot_left = (i == 0) || rot_left;
      break;
    case CTD_RCOAX_WITH_PREV:
      assert(branch > 0 && p[branch] + 1 < int(r.size()));
      // (   ).(   ). or (.(   ).   )
      energy = em.MismatchCoaxial(enb, r[p[branch] + 1], r[branch - 1], stb);
      // Need to do rotations
      branch_ctds.emplace_back(CTD_RCOAX_WITH_NEXT, energy);
      rot_left = (i == 0) || rot_left;
      break;
    case CTD_FCOAX_WITH_PREV:
      // (   )(   ) or ((   )   ) or (   (   ))
      energy = em.stack[r[p[prev_branch]]][stb][enb][r[prev_branch]];
      branch_ctds.emplace_back(CTD_FCOAX_WITH_NEXT, energy);
      rot_left = (i == 0) || rot_left;
      break;
    case CTD_LCOAX_WITH_NEXT:
    case CTD_RCOAX_WITH_NEXT:
    case CTD_FCOAX_WITH_NEXT:
      // All these cases will be handled in the next branch (PREV).
      continue;
    default: verify_expr(false, "bug");  // Should never happen
    }
    branch_ctds.emplace_back(computed.base_ctds[branch], energy);
    total_energy += energy;
  }
  if (rot_left) {
    branch_ctds.push_back(branch_ctds.front());
    branch_ctds.pop_front();
  }
  return total_energy;
}
}  // namespace internal
}  // namespace energy
}  // namespace memerna
