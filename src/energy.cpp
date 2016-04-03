#include <cstdio>
#include <cmath>
#include "energy.h"
#include "globals.h"

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
    printf("Hairpin %d %d applying AUGU penalty\n", st, en);
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
    printf("Applying special c bulge: %d\n", bulge_special_c);
    energy += bulge_special_c;
  }
#if BULGE_LOOP_STATES
  // TODO. This is not implemented and potentially experimentally not super solid.
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

energy_t InternalLoopEnergy(int ost, int oen, int ist, int ien) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 1 && botlen == 1)
    return internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[oen + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[oen + 2]][r[oen + 1]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return internal_1x2[r[oen]][r[oen + 1]][r[oen + 2]][r[ien]][r[ist]][r[ost + 1]][r[ost]];
  if (toplen == 2 && botlen == 2)
    return internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[oen + 2]][r[oen + 1]][r[oen]];
  return 0;
}

energy_t ComputeEnergyInternal(int st, int en) {
  assert(en >= st - 1);
  energy_t energy = 0;
  printf("%d %d\n", st, en);
  // Watson-Crick helix.
  if (p[st + 1] == en - 1 && p[st] == en) {
    printf("Adding Watson-Crick stacking energy %d %d %d %d: %d\n",
           st, st + 1, en - 1, en, stacking_e[r[st]][r[st + 1]][r[en - 1]][r[en]]);
    energy += stacking_e[r[st]][r[st + 1]][r[en - 1]][r[en]];
  }
  // Look for loops inside.
  int numloops = 0;
  int first_st = -1, first_en = -1;
  for (int i = st; i <= en; ++i) {
    int pair = p[i];
    // This only happens when we're called at 0, sz - 1 and it's something like "(...)."
    if ((i == st && pair == en) || (i == en && pair == st)) continue;
    assert(pair < en);
    if (pair != -1) {
      printf("Loop at %d %d\n", i, pair);
      if (first_st == -1) {
        first_st = i;
        first_en = pair;
      }
      energy += ComputeEnergyInternal(i, pair);

      i = pair;
      numloops++;
    }
  }

  // Hairpin loop.
  if (numloops == 0) {
    int numunpaired = en - st - 1;
    assert(numunpaired >= 3);
    printf("Found hairpin loop %d %d\n", st, en);
    energy += HairpinEnergy(st, en);
  } else if (numloops == 1 && first_st - st + en - first_en > 2) {
    // Bulge loop or internal loop.
    if (first_st - st == 1 || en - first_en == 1) {
      printf("Found bulge loop %d %d\n", first_st, first_en);
      energy += BulgeEnergy(st, en, first_st, first_en);
    } else {
      printf("Found internal loop %d %d\n", first_st, first_en);
    }
  }
  printf("Found %d loops\n", numloops);
  return energy;
}

// TODO: This does a copy.
energy_t ComputeEnergy(const folded_rna_t& frna) {
  r = frna.r;
  p = frna.p;
  energy_t firstAu = 0;
  // Add top level AUGU penalties.
  for (int i = 0; i < static_cast<int>(r.size()); ++i) {
    int pair = p[i];
    if (pair != -1) {
      if (IsUnorderedOf(r[i], r[pair], G_b | A_b, U_b))
        firstAu += AUGU_PENALTY;
      i = pair;
    }
  }
  return ComputeEnergyInternal(0, (int) r.size() - 1) + firstAu;
}

}
}
