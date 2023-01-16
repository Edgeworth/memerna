// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_COMMON_T04LIKE_BRANCH_H_
#define COMPUTE_ENERGY_COMMON_T04LIKE_BRANCH_H_

#include <algorithm>
#include <deque>
#include <tuple>
#include <utility>
#include <vector>

#include "compute/energy/common/branch.h"
#include "compute/energy/t04/model.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::erg {

// Computes the optimal arrangement of coaxial stackings, terminal mismatches, and dangles (CTD).
// This DP needs to be run four times. The series of branches is actually cyclic, and there are two
// types of interactions that can occur between branches. Interactions between the branches
// themselves, and interactions between them and the bases next to them. N.B. interact is not
// used in the chemistry sense here. The DP needs to be run with the outer branch at the start
// and the end, and also with the unpaired base on the left of the first branch (if it exists)
// not used and potentially used. Running with the outer branch at both the start and
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
    Energy macro_upd_value_ = (value);                                                \
    if (cache[cur_used][cur_idx] + macro_upd_value_ < cache[used][idx]) {             \
      cache[used][idx] = cache[cur_used][cur_idx] + macro_upd_value_;                 \
      back[used][idx] = std::make_tuple(cur_used, cur_idx, macro_upd_value_, reason); \
    }                                                                                 \
  } while (0)

template <typename T>
requires std::is_base_of_v<T04ModelMixin, T> Energy ComputeOptimalCtds(const T& em,
    const Primary& r, const Secondary& s, const std::deque<int>& branches, bool use_first_lu,
    BranchCtd* branch_ctd) {
  int N = static_cast<int>(branches.size());
  int RSZ = static_cast<int>(r.size());
  assert(branch_ctd->empty());
  // Could be on the exterior loop with a branch (0, N - 1).
  if (N < 1) return ZERO_E;

  // cache[used][i]
  std::vector<Energy> cache[2] = {
      std::vector<Energy>(N + 1, MAX_E), std::vector<Energy>(N + 1, MAX_E)};
  std::vector<std::tuple<bool, int, Energy, Ctd>> back[2] = {
      std::vector<std::tuple<bool, int, Energy, Ctd>>(
          N + 1, std::make_tuple(false, -1, ZERO_E, CTD_NA)),
      std::vector<std::tuple<bool, int, Energy, Ctd>>(
          N + 1, std::make_tuple(false, -1, ZERO_E, CTD_NA))};

  cache[0][0] = cache[1][0] = ZERO_E;
  int first_lui = branches[0] - 1;
  int last_rui = s[branches[N - 1]] + 1;

  // Precompute data about the unpaired bases.
  std::vector<int> li(N);
  std::vector<int> ri(N);
  std::vector<int> lui(N);
  std::vector<int> rui(N);
  std::vector<bool> lu_exists(N);
  std::vector<bool> lu_usable(N);
  std::vector<bool> ru_exists(N);
  std::vector<bool> ru_usable(N);
  std::vector<bool> ru_shared(N);
  for (int i = 0; i < N; ++i) {
    li[i] = branches[i];
    ri[i] = s[branches[i]];
    assert(ri[i] != -1);
    lui[i] = li[i] - 1;
    rui[i] = ri[i] + 1;
    // If |use_first_lu|, then if the left unpaired base is the same as the last branch's right
    // unpaired base, then we can't use it (as it could be used at the end by a terminal mismatch,
    // dangle, right facing coaxial stack, etc). This is because the loop is cyclic.
    lu_exists[i] = lui[i] >= 0 && lui[i] < RSZ && s[lui[i]] == -1;
    lu_usable[i] = lu_exists[i] && (lui[i] != last_rui || use_first_lu);
    ru_exists[i] = rui[i] >= 0 && rui[i] < RSZ && s[rui[i]] == -1;
    ru_usable[i] = ru_exists[i] && (rui[i] != first_lui || !use_first_lu);
    ru_shared[i] = ru_exists[i] && rui[i] < RSZ - 1 && s[rui[i] + 1] != -1;
  }

  for (int i = 0; i < N; ++i) {
    Base lb = r[li[i]];
    Base rb = r[ri[i]];
    Base lub = -1;
    Base rub = -1;
    if (lu_exists[i]) lub = r[lui[i]];
    if (ru_exists[i]) rub = r[rui[i]];

    // Flush coaxial stacking. Requires that ru not exist (i.e. adjacent branches) and this not be
    // the last branch.
    if (!ru_exists[i] && i != N - 1) {
      Energy coax = em.stack[rb][r[li[i + 1]]][r[ri[i + 1]]][lb];
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
        Energy left_coax = em.MismatchCoaxial(rb, rub, lub, lb);
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
        Energy right_coax = em.MismatchCoaxial(r[ri[i + 1]], r[rui[i + 1]], rub, r[li[i + 1]]);

        UPDATE_CACHE(ru_shared[i + 1], i + 2, 0, i, right_coax, CTD_RC_WITH_NEXT);
        if (lu_exists[i]) {
          // In the case that lu doesn't exist but it is "used" it means this branch was consumed by
          // a coaxial interaction so don't use it.
          UPDATE_CACHE(ru_shared[i + 1], i + 2, 1, i, right_coax, CTD_RC_WITH_NEXT);
        }
      }
    }

    if (lu_usable[i]) {
      // 5' dangle.
      UPDATE_CACHE(0, i + 1, 0, i, em.dangle5[rb][lub][lb], CTD_5_DANGLE);
    }

    // Have the option of doing nothing.
    UPDATE_CACHE(0, i + 1, 0, i, ZERO_E, CTD_UNUSED);
    UPDATE_CACHE(0, i + 1, 1, i, ZERO_E, CTD_UNUSED);
  }

  std::tuple<bool, int, Energy, Ctd> state{false, N, ZERO_E, CTD_NA};
  if (cache[1][N] < cache[0][N]) state = std::make_tuple(true, N, ZERO_E, CTD_NA);
  // First state contains no real info, so go ahead one.
  state = back[std::get<0>(state)][std::get<1>(state)];
  assert(branch_ctd->empty());
  while (true) {
    bool used = false;
    int idx = 0;
    Energy energy = ZERO_E;
    Ctd reason = CTD_NA;
    std::tie(used, idx, energy, reason) = std::move(state);
    if (idx == -1) break;
    // We can skip a branch on coaxial stacking interactions, so make sure to insert the ctd energy
    // for the branch. We build branch_ctd backwards, so we need to insert this later branch first.
    // To get the PREV versions, just add one.
    if (reason == CTD_LCOAX_WITH_NEXT || reason == CTD_RC_WITH_NEXT ||
        reason == CTD_FCOAX_WITH_NEXT)
      branch_ctd->emplace_back(Ctd(reason + 1), energy);
    branch_ctd->emplace_back(reason, energy);
    state = back[used][idx];
  }
  std::reverse(branch_ctd->begin(), branch_ctd->end());

  return std::min(cache[0][N], cache[1][N]);
}

#undef UPDATE_CACHE

// Reads the per-base ctd representation from |ctd| for |branches| branches and
// writes it in branch representation to |branch_ctd|.
template <typename T>
requires std::is_base_of_v<T04ModelMixin, T> Energy AddBaseCtdsToBranchCtds(const T& em,
    const Primary& r, const Secondary& s, const Ctds& ctd, const std::deque<int>& branches,
    BranchCtd* branch_ctd) {
  assert(branch_ctd->empty());
  Energy total_energy = ZERO_E;
  // If we have an outer loop in |branches|, it is possible the first could refer to PREV, or the
  // last, to NEXT. In this case, we need to fix the branch_ctd so that the corresponding branch
  // ctd is on the right side. e.g. if the first element refers to PREV, we would put something
  // before it, but that actually needs to be at the end.
  bool rot_left = false;
  for (int i = 0; i < static_cast<int>(branches.size()); ++i) {
    const int branch = branches[i];
    const int prev_branch = i > 0 ? branches[i - 1] : branches.back();
    Energy energy = ZERO_E;
    const auto stb = r[branch];
    const auto enb = r[s[branch]];
    switch (ctd[branch]) {
    case CTD_UNUSED: break;
    case CTD_3_DANGLE:
      assert(s[branch] + 1 < static_cast<int>(r.size()));
      energy = em.dangle3[enb][r[s[branch] + 1]][stb];
      break;
    case CTD_5_DANGLE:
      assert(branch > 0);
      energy = em.dangle5[enb][r[branch - 1]][stb];
      break;
    case CTD_MISMATCH:
      assert(s[branch] + 1 < static_cast<int>(r.size()) && branch > 0);
      energy = em.terminal[enb][r[s[branch] + 1]][r[branch - 1]][stb];
      break;
    case CTD_LCOAX_WITH_PREV:
      // .(   ).(   )
      assert(s[prev_branch] + 1 < static_cast<int>(r.size()) && prev_branch - 1 >= 0);
      energy = em.MismatchCoaxial(
          r[s[prev_branch]], r[s[prev_branch] + 1], r[prev_branch - 1], r[prev_branch]);
      branch_ctd->emplace_back(CTD_LCOAX_WITH_NEXT, energy);
      rot_left = (i == 0) || rot_left;
      break;
    case CTD_RC_WITH_PREV:
      assert(branch > 0 && s[branch] + 1 < static_cast<int>(r.size()));
      // (   ).(   ). or (.(   ).   )
      energy = em.MismatchCoaxial(enb, r[s[branch] + 1], r[branch - 1], stb);
      // Need to do rotations
      branch_ctd->emplace_back(CTD_RC_WITH_NEXT, energy);
      rot_left = (i == 0) || rot_left;
      break;
    case CTD_FCOAX_WITH_PREV:
      // (   )(   ) or ((   )   ) or (   (   ))
      energy = em.stack[r[s[prev_branch]]][stb][enb][r[prev_branch]];
      branch_ctd->emplace_back(CTD_FCOAX_WITH_NEXT, energy);
      rot_left = (i == 0) || rot_left;
      break;
    case CTD_LCOAX_WITH_NEXT:
    case CTD_RC_WITH_NEXT:
    case CTD_FCOAX_WITH_NEXT:
      // All these cases will be handled in the next branch (PREV).
      continue;
    default: bug();  // Should never happen
    }
    branch_ctd->emplace_back(ctd[branch], energy);
    total_energy += energy;
  }
  if (rot_left) {
    branch_ctd->push_back(branch_ctd->front());
    branch_ctd->pop_front();
  }
  return total_energy;
}

}  // namespace mrna::erg

#endif  // COMPUTE_ENERGY_COMMON_T04LIKE_BRANCH_H_
