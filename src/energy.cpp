#include <cstdio>
#include <cmath>
#include "energy.h"
#include "globals.h"

#if ENERGY_LOG
#define ELOG(msg, ...) fprintf(stderr, "L%.4d: %s" msg, __LINE__, __func__, __VA_ARGS__)
#else
#define ELOG(msg, ...)
#endif

namespace memerna {
namespace energy {

energy_t HairpinInitiation(int n) {
  assert(n >= 3);
  if (n < INITIATION_CACHE_SZ) return hairpin_init[n];
  static_assert(INITIATION_CACHE_SZ > 9, "Need initiation values for up to 9.");
  // Formula: G_init(9) + 1.75 * R * T * ln(n / 9).
  return energy_t(round(hairpin_init[9] + 10.0 * 1.75 * R * T * log(n / 9.0)));
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
//   If the mismatch is UU or GA (in that order), additional bonus (TODO: bake this into the mismatch energies for hairpin loops)
//   If the mismatch is GG, additional bonus.
//   If the pair st, en is GU (not UG), a bonus if st - 1 and st - 2 are both Gs.
//   A penalty if all the bases inside are C: A * length + B (A, B specified as part of the energy model).
energy_t HairpinEnergy(int st, int en) {
  assert(st < en);
  std::string s;
  for (int i = st; i <= en; ++i)
    s += BaseToChar(r[i]);
  auto iter = hairpin_e.find(s);
  if (iter != hairpin_e.end())
    return iter->second;

  // Subtract two for the initiating base pair.
  int length = en - st - 1;
  if (length < 3) return MAX_E;  // Disallowed by T04.
  energy_t energy = HairpinInitiation(length);
  // Apply AU penalty if necessary (N.B. not for special hairpin sequences).
  if (IsUnorderedOf(r[st], r[en], G_b | A_b, U_b)) {
    ELOG("(%d, %d): applying AUGU penalty\n", st, en);
    energy += AUGU_PENALTY;
  }

  // T04 says hairpin loops with all C bases inside them are treated specially.
  bool all_c = true;
  for (int i = st + 1; i <= en - 1; ++i) {
    if (r[i] != C) all_c = false;
  }

  if (length == 3)
    return energy + (all_c ? hairpin_c3_loop : 0);
  base_t left = r[st + 1], right = r[en - 1];
  energy += terminal_e[r[st]][left][right][r[en]];
  if (IsPairOf(left, right, U_b, U_b) || IsPairOf(left, right, G_b, A_b))
    energy += hairpin_uu_ga_first_mismatch;
  if (IsPairOf(left, right, G_b, G_b))
    energy += hairpin_gg_first_mismatch;
  if (all_c)
    energy += hairpin_all_c_a * length + hairpin_all_c_b;

  assert(st >= 2);
  if (IsPairOf(r[st], r[en], G_b, U_b) && r[st - 1] == G && r[st - 2] == G)
    energy += hairpin_special_gu_closure;
  return energy;
}

energy_t BulgeInitiation(int n) {
  assert(n >= 1);
  if (n < INITIATION_CACHE_SZ) return bulge_init[n];
  static_assert(INITIATION_CACHE_SZ > 6, "Need initiation values for up to 6.");
  // Formula: G_init(6) + 1.75 * R * T * ln(n / 6).
  return energy_t(round(bulge_init[6] + 10.0 * 1.75 * R * T * log(n / 6.0)));
}

// Indices are inclusive.
// Rules for bulge loop energy:
// 1. If length (number of unpaired bases in the bulge loop) is more than 1, just use initiation energy.
// 2. Otherwise:
//    Don't apply AU/GU penalties -- since the helix is able to continue (unlike for lengths > 1).
//    Since the helix continues, also apply stacking energies for Watson-Crick helices.
//    If the unpaired base is a C, and is next to another C (at pos - 1 or pos + 1), add special C bulge bonus.
energy_t BulgeEnergy(int ost, int oen, int ist, int ien) {
  assert(ist > ost && ien < oen && (oen - ien == 1 || ist - ost == 1) && (oen - ien >= 2 || ist - ost >= 2));
  int length = std::max(ist - ost, oen - ien) - 1;
  energy_t energy = BulgeInitiation(length);
  if (length > 1) {
    // Bulges of length > 1 are considered separate helices and get AU/GU penalties.
    if (IsUnorderedOf(r[ost], r[oen], G_b | A_b, U_b))
      energy += AUGU_PENALTY;
    if (IsUnorderedOf(r[ist], r[ien], G_b | A_b, U_b))
      energy += AUGU_PENALTY;
    return energy;
  }
  // Stacking energy.
  energy += stacking_e[r[ost]][r[ist]][r[ien]][r[ist]];
  int unpaired = ost + 1;
  if (ost + 1 == ist) unpaired = ien + 1;
  // Special C bulge.
  if (r[unpaired] == C && (r[unpaired - 1] == C || r[unpaired] + 1 == C)) {
    ELOG("(%d, %d, %d, %d): applying special c bulge: %d\n", ost, oen, ist, ien, bulge_special_c);
    energy += bulge_special_c;
  }
#if BULGE_LOOP_STATES
  // TODO. This is not implemented and potentially experimentally not super solid?
  // Count up the number of contiguous same bases next to the size 1 bulge loop base.
#endif
  return energy;
}

energy_t InternalLoopInitiation(int n) {
  assert(n >= 4);
  if (n < INITIATION_CACHE_SZ) return internal_init[n];
  static_assert(INITIATION_CACHE_SZ > 6, "Need initiation values for up to 6.");
  // Formula: G_init(6) + 1.08 * ln(n / 6).
  return energy_t(round(internal_init[6] + 10.0 * 1.08 * log(n / 6.0)));
}

// Indices are inclusive.
// Rules for internal loops:
// 1. If it is 1x1, 1x2, 2x1, or 2x2, then check a lookup table.
// 2. If not, return G_init plus:
// 2.1 Internal loop specific AU/GU penalty for each AU/GU end.
// 2.2 A constant times the absolute difference between the number of unpaired bases on each side of the loop.
// 2.3 If the loop is 2x3 or 3x2, look up special mismatch parameters. We just store the values for 2x3, and then
//   rotate the rna by 180 degrees to look it up for 3x2.
energy_t InternalLoopEnergy(int ost, int oen, int ist, int ien) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  ELOG("(%d, %d, %d, %d): Computing %dx%d internal loop\n", ost, oen, ist, ien, toplen, botlen);
  if (toplen == 1 && botlen == 1) {
    ELOG("(%d, %d, %d, %d): 1x1 internal loop\n", ost, oen, ist, ien);
    return internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  }
  if (toplen == 1 && botlen == 2) {
    ELOG("(%d, %d, %d, %d): 1x2 internal loop\n", ost, oen, ist, ien);
    return internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  }
  if (toplen == 2 && botlen == 1) {
    ELOG("(%d, %d, %d, %d): 2x1 internal loop\n", ost, oen, ist, ien);
    return internal_1x2[r[ien]][r[ien + 1]][r[ien + 2]][r[oen]][r[ost]][r[ost + 1]][r[ist]];
  }
  if (toplen == 2 && botlen == 2) {
    ELOG("(%d, %d, %d, %d): 2x2 internal loop\n", ost, oen, ist, ien);
    return internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  }

  energy_t energy = InternalLoopInitiation(toplen + botlen);
  ELOG("(%d, %d, %d, %d): Internal loop initiation: %d\n", ost, oen, ist, ien, energy);
  // Special AU/GU penalties.
  if (IsUnorderedOf(r[ost], r[oen], G_b | A_b, U_b)) {
    ELOG("(%d, %d, %d, %d): Adding outer internal loop AUGU penalty\n", ost, oen, ist, ien);
    energy += internal_augu_penalty;
  }
  if (IsUnorderedOf(r[ist], r[ien], G_b | A_b, U_b)) {
    ELOG("(%d, %d, %d, %d): Adding inner internal loop AUGU penalty\n", ost, oen, ist, ien);
    energy += internal_augu_penalty;
  }
  // Asymmetry term.
  energy += std::abs(toplen - botlen) * internal_asym;

  // Special mismatch parameters.
  // To flip an RNA, we flip it vertically and horizontally (180 degree rotation).
  // It turns out that the accesses for 2x3 and 3x2 both touch the same array locations, just which side is
  // on the left / right is flipped.
  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2)) {
    ELOG("(%d, %d, %d, %d): Adding 2x3 internal loop mismatch parameters.\n", ost, oen, ist, ien);
    energy += internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]];
    energy += internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
  } else if (toplen != 1 && botlen != 1) {
    ELOG("(%d, %d, %d, %d): Adding other size internal loop mismatch parameters.\n", ost, oen, ist, ien);
    energy += internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]];
    energy += internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
  }

  return energy;
}

// 1999 rules.
energy_t MultiloopT99Initiation(int num_unpaired, int num_branches) {
  if (num_unpaired > 6)
    return energy_t(
        multiloop_t99_a + 6 * multiloop_t99_b + 10.0 * 1.1 * log(num_unpaired / 6.0) + multiloop_t99_c * num_branches);
  return energy_t(multiloop_t99_a + multiloop_t99_b * num_unpaired + multiloop_t99_c * num_branches);
}

energy_t MultiloopHackInitiation(int num_branches) {
  return multiloop_hack_a + num_branches * multiloop_hack_b;
}

// Computes the optimal arrangement of coaxial stackings, terminal mismatches, and dangles (CTD).
// This DP needs to be run four times. The series of branches is actually cyclic, and there are two types of interactions
// that can occur between branches. Interactions between the branches themselves, and interactions between them and the
// bases next to them. N.B. interact is not used in the chemistry sense here.
// The DP needs to be run with the outer branch at the start and the end, and also with the unpaired base on the left
// // of the first branch (if it exists) not used and potentially used. Running with the outer branch at both the start and
// the end allows coaxial stacking interactions between the outer branch and both adjacent branches. Running with having the
// left unpaired base both not used lets the last branch potentially consume it (wrapping around).
//
// The state holds which branches we're up to, |i|, and whether the unpaired base on the left was consumed (if it exists).
//
// Note that we can't form two dangles attached to the same stem; that's a terminal mismatch.

// Rules for mismatch mediated coaxial stacking:
//    1. A terminal mismatch is formed around the branch being straddled.
//    2. An arbitrary bonus is added.
//    2. An arbitrary bonus is added if the mismatch is Watson-Crick or GU.
inline energy_t MismatchMediatedCoaxialEnergy(
    base_t fiveTop, base_t mismatch_top, base_t mismatch_bot, base_t threeBot) {
  energy_t coax = terminal_e[fiveTop][mismatch_top][mismatch_bot][threeBot] + coax_mismatch_non_contiguous;
  if (IsUnorderedOf(mismatch_top, mismatch_bot, G_b, C_b) ||
      IsUnorderedOf(mismatch_top, mismatch_bot, A_b, U_b))
    coax += coax_mismatch_wc_bonus;
  if (IsUnorderedOf(mismatch_top, mismatch_bot, G_b, U_b))
    coax += coax_mismatch_gu_bonus;
  return coax;
}

// unpaired base will never have been consumed by a 5' dangle.
#define UPDATE_CACHE(used, idx, cur_used, cur_idx, value, reason) \
  do { if (cache[cur_used][cur_idx] + value < cache[used][idx]) { \
    cache[used][idx] = cache[cur_used][cur_idx] + value; \
    back[used][idx] = std::make_tuple(cur_used, cur_idx, \
      std::string(reason) + " energy: " + std::to_string(value) + " total: " + std::to_string(cache[used][idx])); \
  } } while (0)

energy_t ComputeOptimalCtd(const std::deque<int>& branches, int outer_idx, bool use_first_lu) {
  // cache[used][i]
  int N = branches.size();
  int R = r.size();
  assert(outer_idx == 0 || outer_idx == N - 1 || outer_idx == -1);

  std::vector<int> cache[2] = {
      std::vector<int>(size_t(N + 1), MAX_E),
      std::vector<int>(size_t(N + 1), MAX_E)
  };
  std::vector<std::tuple<bool, int, std::string>> back[2] = {
      std::vector<std::tuple<bool, int, std::string>>(size_t(N + 1), std::make_tuple(false, -1, "")),
      std::vector<std::tuple<bool, int, std::string>>(size_t(N + 1), std::make_tuple(false, -1, ""))
  };

  cache[0][0] = cache[1][0] = 0;
  // These values are for outer_idx == N - 1, i.e. last.
  int first_lui = -1, last_rui = -1;
  // These are for outer_idx == 0.
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
    // If |use_first_lu|, then if the left unpaired base is the same as the last branches right unpaired base,
    // then we can't use it (as it could be used at the end by a terminal mismatch, dangle, right facing coaxial stack,
    // etc). This is because the loop is cyclic.
    lu_exists[i] = lui[i] >= 0 && lui[i] < R && p[lui[i]] == -1;
    lu_usable[i] = lu_exists[i] && (lui[i] != last_rui || use_first_lu);
    ru_exists[i] = rui[i] >= 0 && rui[i] < R && p[rui[i]] == -1;
    ru_usable[i] = ru_exists[i] && (rui[i] != first_lui || !use_first_lu);
    ru_shared[i] = ru_exists[i] && rui[i] < R - 1 && p[rui[i] + 1] != -1;
  }

  for (int i = 0; i < N; ++i) {
    base_t lb = r[li[i]], rb = r[ri[i]], lub = -1, rub = -1;
    if (lu_exists[i]) lub = r[lui[i]];
    if (ru_exists[i]) rub = r[rui[i]];

    // Flush coaxial stacking. Requires that ru not exist (i.e. adjacent branches) and this not be the last branch.
    if (!ru_exists[i] && i != N - 1) {
      energy_t coax = stacking_e[rb][r[li[i + 1]]][r[ri[i + 1]]][lb];
      // When the next branch is consumed by this coaxial stack, it can no longer interact with anything, so
      // just skip to i + 2.
      UPDATE_CACHE(0, i + 2, 0, i, coax, "flush coaxial; no lu;");
      if (lu_exists[i]) {
        // If lu exists, and it was used, then it's fine to coaxially stack. If |used| were true but lu didn't exist then
        // we couldn't coaxially stack as the current branch would already have been involved in one, though.
        UPDATE_CACHE(0, i + 2, 1, i, coax, "flush coaxial; lu used;");
      }
    }

    if (lu_usable[i] && ru_usable[i]) {
      // Terminal mismatch, requires lu_exists, ru_exists, and that we didn't use left.
      // Consumes ru, so if it was shared, use it.
      UPDATE_CACHE(ru_shared[i], i + 1, 0, i, terminal_e[rb][rub][lub][lb], "terminal mismatch;");

      // Mismatch mediated coaxial stacking, left facing (uses the branch we're currently at).
      // Requires lu_usable, ru_usable, ru_shared, and left not used. Consumes ru.
      // Skip to the branch after next since the next branch can't be involved in any more interactions anyway:
      // its left pair is consumed, and its right pair can't dangle towards it.
      if (ru_shared[i] && i != N - 1) {
        energy_t left_coax = MismatchMediatedCoaxialEnergy(rb, rub, lub, lb);
        UPDATE_CACHE(0, i + 2, 0, i, left_coax, "left coaxial;");
      }
    }

    // Right consuming cases.
    if (ru_usable[i]) {
      // Right dangle (3').
      // Only requires ru_exists so handle where left is both used and not used.
      UPDATE_CACHE(ru_shared[i], i + 1, 0, i, dangle3_e[rb][rub][lb], "lu unused; ru 3' dangle;");
      UPDATE_CACHE(ru_shared[i], i + 1, 1, i, dangle3_e[rb][rub][lb], "lu used; ru 3' dangle;");

      // Mismatch mediated coaxial stacking, right facing (uses the next branch).
      // Requires ru_exists, ru_shared. Consumes ru and rru.
      if (ru_shared[i] && i != N - 1 && ru_usable[i + 1]) {
        energy_t right_coax = MismatchMediatedCoaxialEnergy(r[ri[i + 1]], r[rui[i + 1]], rub, r[li[i + 1]]);
        UPDATE_CACHE(1, i + 2, 0, i, right_coax, "right coaxial; no lu;");
        if (lu_exists[i]) {
          // In the case that lu doesn't exist but it is "used" it means this branch was consumed by a coaxial interaction
          // so don't use it.
          UPDATE_CACHE(1, i + 2, 1, i, right_coax, "right coaxial; lu used;");
        }
      }
    }

    // 5' dangle.
    UPDATE_CACHE(0, i + 1, 0, i, dangle5_e[rb][lub][lb], "lu 5' dangle;");

    // Have the option of doing nothing.
    UPDATE_CACHE(0, i + 1, 0, i, 0, "no interaction;");
    UPDATE_CACHE(0, i + 1, 1, i, 0, "no interaction; lu used;");
  }

  std::tuple<bool, int, std::string> state{false, N, ""};
  if (cache[1][N] < cache[0][N])
    state = std::make_tuple(true, N, "");
  while (1) {
    bool used;
    int idx;
    std::string reason;
    std::tie(used, idx, reason) = std::move(state);
    if (idx == -1) break;
    ELOG("(%d, %d): %d %d %s\n", outer_idx, int(use_first_lu), int(used), idx, reason.c_str());
    state = back[used][idx];
  }

  return std::min(cache[0][N], cache[1][N]);
}

#undef UPDATE_CACHE

energy_t MultiloopEnergy(int st, int en, std::deque<int>& branches) {
  bool exterior_loop = st == 0 && en == int(r.size() - 1) && p[st] != en;
  energy_t energy = 0;

  // Add AUGU penalties.
  int num_unpaired = 0;
  for (auto branch_st : branches) {
    num_unpaired += p[branch_st] - branch_st + 1;

    if (IsUnorderedOf(r[branch_st], r[p[branch_st]], G_b | A_b, U_b)) {
      ELOG("(%d, %d): Applying opening AUGU penalty at %d %d\n", st, en, branch_st, p[branch_st - st]);
      energy += AUGU_PENALTY;
    }
  }
  num_unpaired = en - st - 1 - num_unpaired;

  if (exterior_loop) {
    // No initiation for the exterior loop.
    energy_t a = ComputeOptimalCtd(branches, -1, true);
    ELOG("(%d, %d): %d exterior\n", st, en, a);
    energy += a;
  } else {
    if (IsUnorderedOf(r[st], r[en], G_b | A_b, U_b)) {
      ELOG("(%d, %d): Applying closing AUGU penalty\n", st, en);
      energy += AUGU_PENALTY;
    }
    energy_t initiation = MultiloopT99Initiation(num_unpaired, int(branches.size() + 1));
    ELOG("(%d, %d): Initiation: %d\n", st, en, initiation);
    energy += initiation;
    branches.push_front(st);
    energy_t a = ComputeOptimalCtd(branches, 0, true);
    energy_t b = ComputeOptimalCtd(branches, 0, false);
    branches.pop_front();
    branches.push_back(st);
    energy_t c = ComputeOptimalCtd(branches, int(branches.size() - 1), true);
    energy_t d = ComputeOptimalCtd(branches, int(branches.size() - 1), false);
    ELOG("(%d, %d): %d %d %d %d\n", st, en, a, b, c, d);
    energy += std::min(a, std::min(b, std::min(c, d)));
  }

  return energy;
}

energy_t ComputeEnergyInternal(int st, int en) {
  assert(en > st);
  energy_t energy = 0;
  // Watson-Crick helix.
  if (p[st + 1] == en - 1 && p[st] == en) {
    ELOG(
        "(%d, %d): Adding Watson-Crick stacking energy %d %d %d %d: %d\n",
        st, en, st, st + 1, en - 1, en, stacking_e[r[st]][r[st + 1]][r[en - 1]][r[en]]);
    energy += stacking_e[r[st]][r[st + 1]][r[en - 1]][r[en]];
  }

  // Look for branches inside.
  std::deque<int> branches;
  for (int i = st; i <= en; ++i) {
    int pair = p[i];
    assert(pair <= en);
    // If we're in the exterior loop we still want to consider branches starting or ending at st or en.
    if ((i == st && pair == en) || (i == en && pair == st)) continue;

    if (pair != -1) {
      branches.push_back(i);
      energy_t branch_energy = ComputeEnergyInternal(i, pair);
      ELOG("(%d, %d): Branch at %d %d: %d\n", st, en, i, pair, branch_energy);
      energy += branch_energy;
      i = pair;
    }
  }

  // We're in the exterior loop if we were called with the entire RNA and there's no match on the very ends that takes
  // us out of the exterior loop.
  bool exterior_loop = st == 0 && en == int(r.size() - 1) && p[st] != en;
  if (exterior_loop || branches.size() >= 2) {
    // Found multiloop.
    energy_t multiloop_energy = MultiloopEnergy(st, en, branches);
    ELOG("(%d, %d): Found multiloop: %d\n", st, en, multiloop_energy);
    energy += multiloop_energy;
  } else if (branches.size() == 0) {
    // Hairpin loop.
    int num_unpaired = en - st - 1;
    assert(num_unpaired >= 3);
    energy_t hairpin_energy = HairpinEnergy(st, en);
    ELOG("(%d, %d): Found hairpin loop: %d\n", st, en, hairpin_energy);
    energy += hairpin_energy;
  } else if (branches.size() == 1) {
    int loop_st = branches.front(), loop_en = p[branches.front()];
    if (loop_st - st + en - loop_en > 2) {
      // Bulge loop or internal loop.
      if (loop_st - st == 1 || en - loop_en == 1) {
        energy_t bulge_energy = BulgeEnergy(st, en, loop_st, loop_en);
        ELOG("(%d, %d): Found bulge loop %d %d %d %d: %d\n", st, en, st, en, loop_st, loop_en, bulge_energy);
        energy += bulge_energy;
      } else {
        energy_t internal_energy = InternalLoopEnergy(st, en, loop_st, loop_en);
        ELOG("(%d, %d): Found internal loop %d %d %d %d: %d\n", st, en, st, en, loop_st, loop_en, internal_energy);
        energy += internal_energy;
      }
    }
  }

  return energy;
}

// TODO: This does a copy.
energy_t ComputeEnergy(const folded_rna_t& frna) {
  r = frna.r;
  p = frna.p;
  ELOG("(): Computing energy of %d length RNA.\n", (int) r.size());
  energy_t firstAu = 0;
  if (p[0] == int(r.size() - 1) && IsUnorderedOf(r[0], r[p[0]], G_b | A_b, U_b))
    firstAu = AUGU_PENALTY;
  return ComputeEnergyInternal(0, (int) r.size() - 1) + firstAu;
}

}
}
