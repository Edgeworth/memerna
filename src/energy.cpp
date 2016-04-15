#include <cstdio>
#include <cmath>
#include "energy.h"
#include "globals.h"

#if ENERGY_LOG
#define ELOG(msg, args...) fprintf(stderr, "%s" msg, __func__, args)
#else
#define ELOG(msg, args...)
#endif

namespace memerna {
namespace energy {

energy_t HairpinInitiation(int n) {
  assert(n >= 3);
  if (n < INITIATION_CACHE_SZ) return hairpin_init[n];
  static_assert(INITIATION_CACHE_SZ > 9, "Need initiation values for up to 9.");
  // Formula: G_init(9) + 1.75 * R * T * ln(n / 9).
  return static_cast<energy_t>(round(hairpin_init[9] + 1.75 * R * T * log(n / 9.0)));
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
  return static_cast<energy_t>(round(bulge_init[6] + 1.75 * R * T * log(n / 6.0)));
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
  return static_cast<energy_t>(round(internal_init[6] + 1.08 * log(n / 6.0)));
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
energy_t MultiloopInitiation(int st, int en, const std::vector<int>& loops) {
  int num_unpaired = 0;
  for (auto loop_st : loops) {
    num_unpaired += loop_st - p[loop_st] + 1;
  }
  num_unpaired = en - st - 1 - num_unpaired;
#if T99_LN_FACTOR
  if (num_unpaired > 6)
    return static_cast<energy_t>(multiloop_a + 6 * multiloop_b + 1.1 * log(num_unpaired / 6.0) + multiloop_c * loops.size());
#endif
  return static_cast<energy_t>(multiloop_a + multiloop_b * num_unpaired + multiloop_c * loops.size());
}

// Computes the optimal arrangement of coaxial stackings, terminal mismatches, and dangles (CTD).
energy_t ComputeOptimalCtd(std::vector<int>& loops, int outer_idx) {
  int N = static_cast<int>(loops.size());
  std::vector<int> cache[2] = {std::vector<int>(loops.size(), MAX_E), std::vector<int>(loops.size(), MAX_E)};
  for (int i = 0; i < N - 1; ++i) {
    int l = loops[i], r = p[loops[i]];
    assert(r != -1);
    int lu = l - 1, ru = r + 1;
    if (i == outer_idx) {
      lu = r - 1;
      ru = l + 1;
    }
    bool lu_exists = p[lu] == -1, ru_exists = p[ru] == -1;
    bool lu_shared = lu_exists && lu > 0 && p[lu - 1] != -1;
    bool ru_shared = ru_exists && ru < N - 1 && p[ru + 1] != -1;
    // Left consuming cases. Requires lu_exists and left not used.
    if (lu_exists) {
      // Terminal mismatch. Requires lu_exists, ru_exists, and left not used.
      cache[ru_shared][i + 1] = std::min(cache[ru_shared][i + 1], cache[0][i] + terminal_e[r][ru][lu][l]);
    }

    // Cases where the left was consumed.

    // C
  }
}

energy_t MultiloopEnergy(int st, int en, const std::vector<int>& loops) {
  bool exterior_loop = st == 0 && en == static_cast<int>(r.size() - 1) && p[st] != en;
  energy_t energy = 0;
  // No initiation for the exterior loop.
  if (!exterior_loop)
    energy += MultiloopInitiation(st, en, loops);

  // Add AUGU penalties.
  for (auto loop_st : loops) {
    if (IsUnorderedOf(r[loop_st], r[p[loop_st]], G_b | A_b, U_b)) {
      ELOG("(%d, %d): Applying opening AUGU penalty at %d %d\n", st, en, loop_st, p[loop_st - st]);
      energy += AUGU_PENALTY;
    }
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

  // Look for loops inside.
  std::vector<int> loops;
  for (int i = st; i <= en; ++i) {
    int pair = p[i];
    assert(pair <= en);
    // If we're in the exterior loop we still want to consider loops starting or ending at st or en.
    if ((i == st && pair == en) || (i == en && pair == st)) continue;

    if (pair != -1) {
      ELOG("(%d, %d): Loop at %d %d\n", st, en, i, pair);
      loops.push_back(i);
      energy += ComputeEnergyInternal(i, pair);

      i = pair;
    }
  }

  // We're in the exterior loop if we were called with the entire RNA and there's no match on the very ends that takes
  // us out of the exterior loop.
  bool exterior_loop = st == 0 && en == static_cast<int>(r.size() - 1) && p[st] != en;
  if (exterior_loop || loops.size() >= 2) {
    // Found multiloop.
    energy_t multiloop_energy = MultiloopEnergy(st, en, loops);
    ELOG("(%d, %d): Found multiloop: %d\n", st, en, multiloop_energy);
    energy += multiloop_energy;
  } else if (loops.size() == 0) {
    // Hairpin loop.
    int num_unpaired = en - st - 1;
    assert(num_unpaired >= 3);
    energy_t hairpin_energy = HairpinEnergy(st, en);
    ELOG("(%d, %d): Found hairpin loop: %d\n", st, en, hairpin_energy);
    energy += hairpin_energy;
  } else if (loops.size() == 1) {
    int loop_st = loops.front(), loop_en = p[loops.front()];
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

  ELOG("(%d, %d): Found %d loops\n", st, en, (int) loops.size());
  return energy;
}

// TODO: This does a copy.
energy_t ComputeEnergy(const folded_rna_t& frna) {
  r = frna.r;
  p = frna.p;
  ELOG("(): Computing energy of %d length RNA.\n", (int)r.size());
  energy_t firstAu = 0;
  if (p[0] == static_cast<int>(r.size() - 1) && IsUnorderedOf(r[0], r[p[0]], G_b | A_b, U_b))
    firstAu = AUGU_PENALTY;
  return ComputeEnergyInternal(0, (int) r.size() - 1) + firstAu;
}

}
}
