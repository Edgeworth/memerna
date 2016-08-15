#include <cstdio>
#include <cmath>
#include <memory>
#include "energy.h"
#include "structure.h"
#include "constants.h"

namespace memerna {
namespace energy {

energy_t HairpinInitiation(int n) {
  assert(n >= 3);
  if (n < INITIATION_CACHE_SZ) return hairpin_init[n];
  static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
  // Formula: G_init(9) + 1.75 * R * T * ln(n / 9)  -- we use 30 here though to match RNAstructure.
  return energy_t(round(hairpin_init[30] + 10.0 * 1.75 * constants::R * constants::T * log(n / 30.0)));
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
//   If the mismatch is UU or GA (not AG), additional bonus (TODO: bake this into the mismatch energies for hairpin loops)
//   If the mismatch is GG, additional bonus.
//   If the pair st, en is GU (not UG), a bonus if st - 1 and st - 2 are both Gs, if they exist.
//   A penalty if all the bases inside are C: A * length + B (A, B specified as part of the energy model).
energy_t HairpinEnergy(int st, int en, std::unique_ptr<structure::Structure>* s) {
  assert(st < en);
  if (s) *s = std::make_unique<structure::HairpinLoop>(st, en);

  std::string seq;
  for (int i = st; i <= en; ++i)
    seq += BaseToChar(r[i]);
  auto iter = hairpin_e.find(seq);
  if (iter != hairpin_e.end())
    return iter->second;

  // Subtract two for the initiating base pair.
  int length = en - st - 1;
  if (length < 3) return constants::MAX_E;  // Disallowed by T04.
  energy_t energy = HairpinInitiation(length);
  if (s) (*s)->AddNote("%de - initiation", energy);
  // Apply AU penalty if necessary (N.B. not for special hairpin sequences).
  if (IsAuGu(r[st], r[en])) {
    if (s) (*s)->AddNote("%de - AU/GU penalty", augu_penalty);
    energy += augu_penalty;
  }

  // T04 says hairpin loops with all C bases inside them are treated specially.
  bool all_c = true;
  for (int i = st + 1; i <= en - 1; ++i) {
    if (r[i] != C) all_c = false;
  }

  if (length == 3) {
    if (all_c) {
      if (s) (*s)->AddNote("%de - all C penalty (length 3)", hairpin_c3_loop);
      energy += hairpin_c3_loop;
    }
    return energy;
  }
  base_t left = r[st + 1], right = r[en - 1];
  if (s) (*s)->AddNote("%de - terminal mismatch", terminal_e[r[st]][left][right][r[en]]);
  energy += terminal_e[r[st]][left][right][r[en]];
  if (IsPairOf(left, right, U_b, U_b) || IsPairOf(left, right, G_b, A_b)) {
    if (s) (*s)->AddNote("%de - UU/GA first mismatch", hairpin_uu_ga_first_mismatch);
    energy += hairpin_uu_ga_first_mismatch;
  }
  if (IsPairOf(left, right, G_b, G_b)) {
    if (s) (*s)->AddNote("%de - GG first mismatch", hairpin_gg_first_mismatch);
    energy += hairpin_gg_first_mismatch;
  }
  if (all_c) {
    energy_t all_c_energy = hairpin_all_c_a * length + hairpin_all_c_b;
    if (s) (*s)->AddNote("%de - all C penalty", all_c_energy);
    energy += all_c_energy;
  }

  if (IsPairOf(r[st], r[en], G_b, U_b) && st >= 2 && r[st - 1] == G && r[st - 2] == G) {
    if (s) (*s)->AddNote("%de - special GU closure", hairpin_special_gu_closure);
    energy += hairpin_special_gu_closure;
  }
  return energy;
}

energy_t BulgeInitiation(int n) {
  assert(n >= 1);
  if (n < INITIATION_CACHE_SZ) return bulge_init[n];
  static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
  // Formula: G_init(6) + 1.75 * R * T * ln(n / 6) -- we use 30 here though to match RNAstructure.
  return energy_t(round(bulge_init[30] + 10.0 * 1.75 * constants::R * constants::T * log(n / 30.0)));
}

// Indices are inclusive.
// Rules for bulge loop energy:
// 1. If length (number of unpaired bases in the bulge loop) is more than 1, just use initiation energy.
// 2. Otherwise:
//    Don't apply AU/GU penalties -- since the helix is able to continue (unlike for lengths > 1).
//    Since the helix continues, also apply stacking energies for Watson-Crick helices.
//    If the unpaired base is a C, and is next to another C (at pos - 1 or pos + 1), add special C bulge bonus.
//    Count up the number of contiguous bases next to the size 1 bulge loop base, and compute a bonus from that.
energy_t BulgeEnergy(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s) {
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
      if (s) (*s)->AddNote("%de - outer AU/GU penalty", augu_penalty);
      energy += augu_penalty;
    }
    if (IsAuGu(r[ist], r[ien])) {
      if (s) (*s)->AddNote("%de - inner AU/GU penalty", augu_penalty);
      energy += augu_penalty;
    }
    return energy;
  }
  // Stacking energy.
  energy += stacking_e[r[ost]][r[ist]][r[ien]][r[oen]];
  int unpaired = ost + 1;
  if (ost + 1 == ist) unpaired = ien + 1;
  // Special C bulge.
  if (r[unpaired] == C && (r[unpaired - 1] == C || r[unpaired + 1] == C)) {
    if (s) (*s)->AddNote("%de - special c bulge", bulge_special_c);
    energy += bulge_special_c;
  }
  // TODO: make this faster?
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

energy_t InternalLoopInitiation(int n) {
  assert(n >= 4);
  if (n < INITIATION_CACHE_SZ) return internal_init[n];
  static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
  // Formula: G_init(6) + 1.08 * ln(n / 6) -- we use 30 here though to match RNAstructure.
  return energy_t(round(internal_init[30] + 10.0 * 1.08 * log(n / 30.0)));
}

// Indices are inclusive.
// Rules for internal loops:
// 1. If it is 1x1, 1x2, 2x1, or 2x2, then check a lookup table.
// 2. If not, return G_init plus:
// 2.1 Internal loop specific AU/GU penalty for each AU/GU end.
// 2.2 A constant times the absolute difference between the number of unpaired bases on each side of the loop.
// 2.3 If the loop is 2x3 or 3x2, look up special mismatch parameters. We just store the values for 2x3, and then
//   rotate the rna by 180 degrees to look it up for 3x2.
energy_t InternalLoopEnergy(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (s) {
    *s = std::make_unique<structure::InternalLoop>(ost, oen, ist, ien);
    (*s)->AddNote("%dx%d internal loop", toplen, botlen);
  }
  if (toplen == 1 && botlen == 1)
    return internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];

  energy_t energy = InternalLoopInitiation(toplen + botlen);
  if (s) (*s)->AddNote("%de - initiation", energy);

  // Special AU/GU penalties.
  if (IsAuGu(r[ost], r[oen])) {
    if (s) (*s)->AddNote("%de - outer AU/GU penalty", internal_augu_penalty);
    energy += internal_augu_penalty;
  }
  if (IsAuGu(r[ist], r[ien])) {
    if (s) (*s)->AddNote("%de - inner AU/GU penalty", internal_augu_penalty);
    energy += internal_augu_penalty;
  }
  // Asymmetry term, limit with Ninio maximum asymmetry.
  energy_t asym = std::min(std::abs(toplen - botlen) * internal_asym, constants::NINIO_MAX_ASYM);
  if (s) (*s)->AddNote("%de - asymmetry", asym);
  energy += asym;

  // Special mismatch parameters.
  // To flip an RNA, we flip it vertically and horizontally (180 degree rotation).
  // It turns out that the accesses for 2x3 and 3x2 both touch the same array locations, just which side is
  // on the left / right is flipped.
  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2)) {
    energy_t mismatch =
        internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
            internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("%de - 2x3 mismatch params", mismatch);
    energy += mismatch;
  } else if (toplen != 1 && botlen != 1) {
    energy_t mismatch =
        internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
            internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("%de - other mismatch params", mismatch);
    energy += mismatch;
  }

  return energy;
}

energy_t TwoLoop(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) {
    if (s) *s = std::make_unique<structure::Stacking>(ost, oen);
    return stacking_e[r[ost]][r[ist]][r[ien]][r[oen]];
  }
  if (toplen >= 1 && botlen >= 1)
    return InternalLoopEnergy(ost, oen, ist, ien, s);
  return BulgeEnergy(ost, oen, ist, ien, s);
}

energy_t MultiloopInitiation(int num_branches) {
  return multiloop_hack_a + num_branches * multiloop_hack_b;
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
      back[used][idx] = std::make_tuple(cur_used, cur_idx, \
        std::string(reason) + " energy: " + std::to_string(macro_upd_value_) + \
        " total: " + std::to_string(cache[used][idx])); \
    } \
  } while (0)

energy_t ComputeOptimalCtd(const std::deque<int>& branches, int outer_idx, bool use_first_lu,
    std::unique_ptr<structure::Structure>* s) {
  int N = int(branches.size());
  int RSZ = int(r.size());
  assert(outer_idx == 0 || outer_idx == N - 1 || outer_idx == -1);
  assert(N >= 3 || outer_idx == -1);
  if (N < 1) return 0;

  // cache[used][i]
  std::vector<int> cache[2] = {
      std::vector<int>(size_t(N + 1), constants::MAX_E),
      std::vector<int>(size_t(N + 1), constants::MAX_E)
  };
  std::vector<std::tuple<bool, int, std::string>> back[2] = {
      std::vector<std::tuple<bool, int, std::string>>(size_t(N + 1), std::make_tuple(false, -1, "")),
      std::vector<std::tuple<bool, int, std::string>>(size_t(N + 1), std::make_tuple(false, -1, ""))
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

        UPDATE_CACHE(ru_shared[i + 1], i + 2, 0, i, right_coax, "right coaxial; no lu;");
        if (lu_exists[i]) {
          // In the case that lu doesn't exist but it is "used" it means this branch was consumed by a coaxial interaction
          // so don't use it.
          UPDATE_CACHE(ru_shared[i + 1], i + 2, 1, i, right_coax, "right coaxial; lu used;");
        }
      }
    }

    if (lu_usable[i]) {
      // 5' dangle.
      UPDATE_CACHE(0, i + 1, 0, i, dangle5_e[rb][lub][lb], "lu 5' dangle;");
    }

    // Have the option of doing nothing.
    UPDATE_CACHE(0, i + 1, 0, i, 0, "no interaction;");
    UPDATE_CACHE(0, i + 1, 1, i, 0, "no interaction; lu used;");
  }

  std::tuple<bool, int, std::string> state{false, N, ""};
  if (cache[1][N] < cache[0][N])
    state = std::make_tuple(true, N, "");
  // TODO use ctd representation
  if (s) {
    (*s)->AddNote(
        "%de - CTD: outer_idx: %d, use_first_lu: %d",
        std::min(cache[0][N], cache[1][N]), outer_idx, int(use_first_lu));
  }
  while (1) {
    bool used;
    int idx;
    std::string reason;
    std::tie(used, idx, reason) = std::move(state);
    if (idx == -1) break;
    if (s) (*s)->AddNote("%d %d %s", int(used), idx, reason.c_str());
    state = back[used][idx];
  }

  return std::min(cache[0][N], cache[1][N]);
}

#undef UPDATE_CACHE

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
      if (s) (*s)->AddNote("%de - opening AU/GU penalty at %d %d", augu_penalty, branch_st, p[branch_st]);
      energy += augu_penalty;
    }
  }
  num_unpaired = en - st - 1 - num_unpaired;

  if (exterior_loop) {
    // No initiation for the exterior loop.
    energy += ComputeOptimalCtd(branches, -1, true, s);
    num_unpaired += 2;
  } else {
    if (IsAuGu(r[st], r[en])) {
      if (s) (*s)->AddNote("%de - closing AU/GU penalty at %d %d", augu_penalty, st, en);
      energy += augu_penalty;
    }
    energy_t initiation = MultiloopInitiation(int(branches.size() + 1));
    if (s) (*s)->AddNote("%de - initiation", initiation);
    energy += initiation;
    branches.push_front(st);
    energy_t a = ComputeOptimalCtd(branches, 0, true, s);
    energy_t b = ComputeOptimalCtd(branches, 0, false, s);
    branches.pop_front();
    branches.push_back(st);
    energy_t c = ComputeOptimalCtd(branches, int(branches.size() - 1), true, s);
    energy_t d = ComputeOptimalCtd(branches, int(branches.size() - 1), false, s);
    branches.pop_back();
    if (s) (*s)->AddNote("CTDs: %de %de %de %de", a, b, c, d);
    energy += std::min(a, std::min(b, std::min(c, d)));
  }
  if (s) (*s)->AddNote("Unpaired: %d, Branches: %d", num_unpaired, int(branches.size() + 1));

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
    energy += MultiloopEnergy(st, en, branches, s);
  } else if (branches.size() == 0) {
    // Hairpin loop.
    assert(en - st - 1 >= 3);
    energy += HairpinEnergy(st, en, s);
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

energy_t ComputeEnergy(std::unique_ptr<structure::Structure>* s) {
  energy_t energy = ComputeEnergyInternal(0, (int) r.size() - 1, s);
  if (p[0] == int(r.size() - 1) && IsAuGu(r[0], r[p[0]])) {
    energy += augu_penalty;
    if (s) {
      (*s)->AddNote("%de - top level AU/GU penalty", augu_penalty);
      (*s)->SetSelfEnergy((*s)->GetSelfEnergy() + augu_penalty);  // Gross.
      (*s)->SetTotalEnergy((*s)->GetTotalEnergy() + augu_penalty);  // Gross.
    }
  }
  return energy;
}

}
}
