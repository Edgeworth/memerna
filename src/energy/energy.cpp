#include <cstdio>
#include <cmath>
#include <memory>
#include <algorithm>
#include "energy.h"
#include "structure.h"

namespace memerna {
namespace energy {

namespace {
// Paired array, for easy access.
std::vector<int> p;
std::vector<Ctd> ctds;
}

// Indices are inclusive, include the initiating base pair.
// N.B. This includes an ending AU/GU penalty.
// Rules for hairpin energy:
// 1. Check a lookup table for a hardcoded value.
// 2. Loops with less than 3 bases on the inside are disallowed by T04.
// Case 1 -- 3 bases on the inside:
//   Return G_init + penalty if all bases inside are C.
//   The penalty for all C bases for size 3 hairpin loops is special cased.
// Case 2 -- more than 3 bases on the inside:
//   Return G_init plus the following penalties / bonuses:
//   Terminal mismatch energy for st + 1 and en - 1.
//   If the mismatch is UU or GA (not AG), additional bonus
//   If the mismatch is GG, additional bonus.
//   If the pair st, en is GU (not UG), a bonus if st - 1 and st - 2 are both Gs, if they exist.
//   A penalty if all the bases inside are C: A * length + B (A, B specified as part of the energy model).
energy_t Hairpin(int st, int en, std::unique_ptr<structure::Structure>* s) {
  assert(st < en);
  if (s) *s = std::make_unique<structure::HairpinLoop>(st, en);

  std::string seq;
  for (int i = st; i <= en; ++i)
    seq += BaseToChar(r[i]);
  auto iter = g_hairpin.find(seq);
  if (iter != g_hairpin.end()) {
    if (s) (*s)->AddNote("special hairpin");
    return iter->second;
  }

  // Subtract two for the initiating base pair.
  int length = en - st - 1;
  if (length < 3) return constants::MAX_E;  // Disallowed by T04.
  energy_t energy = HairpinInitiation(length);
  if (s) (*s)->AddNote("%de - initiation", energy);
  // Apply AU penalty if necessary (N.B. not for special hairpin sequences).
  if (IsAuGu(r[st], r[en])) {
    if (s) (*s)->AddNote("%de - AU/GU penalty", g_augu_penalty);
    energy += g_augu_penalty;
  }

  // T04 says hairpin loops with all C bases inside them are treated specially.
  bool all_c = true;
  for (int i = st + 1; i <= en - 1; ++i) {
    if (r[i] != C) all_c = false;
  }

  if (length == 3) {
    if (all_c) {
      if (s) (*s)->AddNote("%de - all C penalty (length 3)", g_hairpin_c3_loop);
      energy += g_hairpin_c3_loop;
    }
    return energy;
  }
  base_t left = r[st + 1], right = r[en - 1];
  if (s) (*s)->AddNote("%de - terminal mismatch", g_terminal[r[st]][left][right][r[en]]);
  energy += g_terminal[r[st]][left][right][r[en]];
  if (IsPairOf(left, right, U_b, U_b) || IsPairOf(left, right, G_b, A_b)) {
    if (s) (*s)->AddNote("%de - UU/GA first mismatch", g_hairpin_uu_ga_first_mismatch);
    energy += g_hairpin_uu_ga_first_mismatch;
  }
  if (IsPairOf(left, right, G_b, G_b)) {
    if (s) (*s)->AddNote("%de - GG first mismatch", g_hairpin_gg_first_mismatch);
    energy += g_hairpin_gg_first_mismatch;
  }
  if (all_c) {
    energy_t all_c_energy = g_hairpin_all_c_a * length + g_hairpin_all_c_b;
    if (s) (*s)->AddNote("%de - all C penalty", all_c_energy);
    energy += all_c_energy;
  }

  if (IsPairOf(r[st], r[en], G_b, U_b) && st >= 2 && r[st - 1] == G && r[st - 2] == G) {
    if (s) (*s)->AddNote("%de - special GU closure", g_hairpin_special_gu_closure);
    energy += g_hairpin_special_gu_closure;
  }
  return energy;
}

// Indices are inclusive.
// Rules for bulge loop energy:
// 1. If length (number of unpaired bases in the bulge loop) is more than 1, just use initiation energy.
// 2. Otherwise:
//    Don't apply AU/GU penalties -- since the helix is able to continue (unlike for lengths > 1).
//    Since the helix continues, also apply stacking energies for Watson-Crick helices.
//    If the unpaired base is a C, and is next to another C (at pos - 1 or pos + 1), add special C bulge bonus.
//    Count up the number of contiguous bases next to the size 1 bulge loop base, and compute a bonus from that.
energy_t Bulge(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s) {
  assert(ist > ost && ien < oen && (oen - ien == 1 || ist - ost == 1) && (oen - ien >= 2 || ist - ost >= 2));
  int length = std::max(ist - ost, oen - ien) - 1;
  energy_t energy = BulgeInitiation(length);

  if (s) {
    *s = std::make_unique<structure::InternalLoop>(ost, oen, ist, ien);
    (*s)->AddNote("Bulge loop len %d", length);
    (*s)->AddNote("%de - initiation", energy);
  }

  if (length > 1) {
    // Bulges of length > 1 are considered separate helices and get AU/GU penalties.
    if (IsAuGu(r[ost], r[oen])) {
      if (s) (*s)->AddNote("%de - outer AU/GU penalty", g_augu_penalty);
      energy += g_augu_penalty;
    }
    if (IsAuGu(r[ist], r[ien])) {
      if (s) (*s)->AddNote("%de - inner AU/GU penalty", g_augu_penalty);
      energy += g_augu_penalty;
    }
    return energy;
  }
  // Stacking energy.
  energy += g_stack[r[ost]][r[ist]][r[ien]][r[oen]];
  int unpaired = ost + 1;
  if (ost + 1 == ist) unpaired = ien + 1;
  // Special C bulge.
  if (r[unpaired] == C && (r[unpaired - 1] == C || r[unpaired + 1] == C)) {
    if (s) (*s)->AddNote("%de - special c bulge", g_bulge_special_c);
    energy += g_bulge_special_c;
  }

  // Count up the number of contiguous same bases next to the size 1 bulge loop base.
  int num_states = 0;
  for (int i = unpaired; i < int(r.size()) && r[i] == r[unpaired]; ++i)
    num_states++;
  for (int i = unpaired - 1; i >= 0 && r[i] == r[unpaired]; --i)
    num_states++;
  energy_t states_bonus = -energy_t(round(10.0 * constants::R * constants::T * log(num_states)));
  if (s) (*s)->AddNote("%de - %d states bonus", states_bonus, num_states);
  energy += states_bonus;

  return energy;
}

// Indices are inclusive.
// Rules for internal loops:
// 1. If it is 1x1, 1x2, 2x1, or 2x2, then check a lookup table.
// 2. If not, return G_init plus:
// 2.1 Internal loop specific AU/GU penalty for each AU/GU end.
// 2.2 A constant times the absolute difference between the number of unpaired bases on each side of the loop.
// 2.3 If the loop is 2x3 or 3x2, look up special mismatch parameters. We just store the values for 2x3, and then
//   rotate the rna by 180 degrees to look it up for 3x2.
energy_t InternalLoop(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (s) {
    *s = std::make_unique<structure::InternalLoop>(ost, oen, ist, ien);
    (*s)->AddNote("%dx%d internal loop", toplen, botlen);
  }
  if (toplen == 1 && botlen == 1)
    return g_internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return g_internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return g_internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return g_internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];

  energy_t energy = InternalLoopInitiation(toplen + botlen);
  if (s) (*s)->AddNote("%de - initiation", energy);

  // Special AU/GU penalties.
  if (IsAuGu(r[ost], r[oen])) {
    if (s) (*s)->AddNote("%de - outer AU/GU penalty", g_internal_augu_penalty);
    energy += g_internal_augu_penalty;
  }
  if (IsAuGu(r[ist], r[ien])) {
    if (s) (*s)->AddNote("%de - inner AU/GU penalty", g_internal_augu_penalty);
    energy += g_internal_augu_penalty;
  }
  // Asymmetry term, limit with Ninio maximum asymmetry.
  energy_t asym = std::min(std::abs(toplen - botlen) * g_internal_asym, constants::NINIO_MAX_ASYM);
  if (s) (*s)->AddNote("%de - asymmetry", asym);
  energy += asym;

  // Special mismatch parameters.
  // To flip an RNA, we flip it vertically and horizontally (180 degree rotation).
  // It turns out that the accesses for 2x3 and 3x2 both touch the same array locations, just which side is
  // on the left / right is flipped.
  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2)) {
    energy_t mismatch =
        g_internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
            g_internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("%de - 2x3 mismatch params", mismatch);
    energy += mismatch;
  } else if (toplen != 1 && botlen != 1) {
    energy_t mismatch =
        g_internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
            g_internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("%de - other mismatch params", mismatch);
    energy += mismatch;
  }

  return energy;
}

energy_t TwoLoop(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) {
    if (s) *s = std::make_unique<structure::Stacking>(ost, oen);
    return g_stack[r[ost]][r[ist]][r[ien]][r[oen]];
  }
  if (toplen >= 1 && botlen >= 1)
    return InternalLoop(ost, oen, ist, ien, s);
  return Bulge(ost, oen, ist, ien, s);
}

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

typedef std::deque<std::pair<Ctd, energy_t>> branch_ctd_t;

// TODO remove outer_idx, do reversal in calling code
energy_t ComputeOptimalCtd(const std::deque<int>& branches, int outer_idx, bool use_first_lu,
    branch_ctd_t* ctd_energies) {
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
        energy_t left_coax = MismatchCoaxial(rb, rub, lub, lb);
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
        energy_t right_coax = MismatchCoaxial(r[ri[i + 1]], r[rui[i + 1]], rub, r[li[i + 1]]);

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
      case CTD_UNUSED: break;
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

energy_t MultiloopEnergy(int st, int en, std::deque<int>& branches, std::unique_ptr<structure::Structure>* s) {
  bool exterior_loop = st == 0 && en == int(r.size() - 1) && p[st] != en;
  energy_t energy = 0;

  if (s) {
    *s = std::make_unique<structure::MultiLoop>(st, en);
    if (exterior_loop) (*s)->AddNote("exterior loop");
  }

  // Add AUGU penalties.
  int num_unpaired = 0;
  for (auto branch_st : branches) {
    num_unpaired += p[branch_st] - branch_st + 1;

    if (IsAuGu(r[branch_st], r[p[branch_st]])) {
      if (s) (*s)->AddNote("%de - opening AU/GU penalty at %d %d", g_augu_penalty, branch_st, p[branch_st]);
      energy += g_augu_penalty;
    }
  }
  num_unpaired = en - st - 1 - num_unpaired + exterior_loop * 2;
  if (s) (*s)->AddNote("Unpaired: %d, Branches: %zu", num_unpaired, branches.size() + 1);

  bool compute_ctd = branches.empty() || (ctds[branches.front()] == CTD_NA);
  branch_ctd_t ctd_energies;
  energy_t ctd_energy = 0;
  if (exterior_loop) {
    // No initiation for the exterior loop.
    if (compute_ctd) {
      ctd_energy = ComputeOptimalCtd(branches, -1, true, &ctd_energies);
      WriteCtdsForBranches(branches, ctd_energies);
    } else {
      ReadCtdsForBranches(branches, &ctd_energies);
    }
  } else {
    if (IsAuGu(r[st], r[en])) {
      if (s) (*s)->AddNote("%de - closing AU/GU penalty at %d %d", g_augu_penalty, st, en);
      energy += g_augu_penalty;
    }
    energy_t initiation = MultiloopInitiation(int(branches.size() + 1));
    if (s) (*s)->AddNote("%de - initiation", initiation);
    energy += initiation;

    if (compute_ctd) {
      branch_ctd_t config_ctds[4] = {};
      std::pair<energy_t, int> config_energies[4] = {};
      branches.push_front(st);
      config_energies[0] = {ComputeOptimalCtd(branches, 0, true, config_ctds), 0};
      config_energies[1] = {ComputeOptimalCtd(branches, 0, false, config_ctds + 1), 1};
      branches.pop_front();
      branches.push_back(st);
      config_energies[2] = {ComputeOptimalCtd(branches, int(branches.size() - 1), true, config_ctds + 2), 2};
      config_energies[3] = {ComputeOptimalCtd(branches, int(branches.size() - 1), false, config_ctds + 3), 3};
      // Shuffle these around so the outer loop is as the start.
      config_ctds[2].push_front(config_ctds[2].back());
      config_ctds[2].pop_back();
      config_ctds[3].push_front(config_ctds[3].back());
      config_ctds[3].pop_back();
      branches.pop_back();
      std::sort(config_energies, config_energies + 4);
      ctd_energies = config_ctds[config_energies[0].second];
      ctd_energy = config_energies[0].first;

      // Assumes the outer loop is the first one.
//      branches.push_front(st);
//      WriteCtdsForBranches(branches, ctd_energies);
//#ifndef NDEBUG
//      branch_ctd_t tmp;
//      assert(ctd_energy == ReadCtdsForBranches(branches, &tmp));
//      assert(ctd_energies == tmp);
//#endif
//      branches.pop_front();
    } else {
      branches.push_front(st);
      ReadCtdsForBranches(branches, &ctd_energies);
      branches.pop_front();
    }
  }
  energy += ctd_energy;

  if (s) {
    (*s)->AddNote("%de - ctd", ctd_energy);
//    int idx = 0;
//    if (!exterior_loop) {
//      (*s)->AddNote("%de - outer loop stacking - %s", ctd_energies[0].second,
//          structure::CtdToName(ctd_energies[0].first));
//      idx++;
//    }
//    for (; idx < ctd_energies.size(); ++i) {
//      (*s)->AddBranch()
//    }
    // TODO note for outer loop stacking
  }

  return energy;
}

energy_t ComputeEnergyInternal(int st, int en, std::unique_ptr<structure::Structure>* s) {
  assert(en >= st);
  energy_t energy = 0;

  // Look for branches inside.
  std::deque<int> branches;
  for (int i = st; i <= en; ++i) {
    int pair = p[i];
    assert(pair <= en);
    if (!(i == st && pair == en) && !(i == en && pair == st) && pair != -1) {
      branches.push_back(i);
      // Skip ahead.
      i = pair;
    }
  }

  // We're in the exterior loop if we were called with the entire RNA and there's no match on the very ends that takes
  // us out of the exterior loop.
  bool exterior_loop = st == 0 && en == int(r.size() - 1) && p[st] != en;
  if (exterior_loop || branches.size() >= 2) {
    // Multiloop.
    energy += MultiloopEnergy(st, en, branches, s);  // TODO add branch occurs below here.
  } else if (branches.size() == 0) {
    // Hairpin loop.
    assert(en - st - 1 >= 3);
    energy += Hairpin(st, en, s);
  } else if (branches.size() == 1) {
    int loop_st = branches.front(), loop_en = p[branches.front()];
    energy += TwoLoop(st, en, loop_st, loop_en, s);
  }

  if (s) (*s)->SetSelfEnergy(energy);
  // Add energy from children.
  for (auto i : branches) {
    int pair = p[i];
    if (s) {
      std::unique_ptr<structure::Structure> structure;
      energy += ComputeEnergyInternal(i, pair, &structure);
      (*s)->AddBranch(std::move(structure));
    } else {
      energy += ComputeEnergyInternal(i, pair, nullptr);
    }
  }
  if (s) (*s)->SetTotalEnergy(energy);

  return energy;
}

computed_t ComputeEnergy(const secondary_t& secondary, std::unique_ptr<structure::Structure>* s) {
  r = secondary.r;
  p = secondary.p;
  ctds.clear();
  ctds.resize(r.size(), CTD_NA);
  assert(r.size() == p.size() && r.size() == ctds.size());
  energy_t energy = ComputeEnergyInternal(0, (int) r.size() - 1, s);
  if (p[0] == int(r.size() - 1) && IsAuGu(r[0], r[p[0]])) {
    energy += g_augu_penalty;
    if (s) {
      (*s)->AddNote("%de - top level AU/GU penalty", g_augu_penalty);
      (*s)->SetSelfEnergy((*s)->GetSelfEnergy() + g_augu_penalty);  // Gross.
      (*s)->SetTotalEnergy((*s)->GetTotalEnergy() + g_augu_penalty);  // Gross.
    }
  }
  return {r, p, ctds, energy};
}

}
}
