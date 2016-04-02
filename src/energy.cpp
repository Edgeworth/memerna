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
// N.B. This doesn't deal with AU/GU penalties.
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
  auto i = hairpin_e.find(s);
  if (i != hairpin_e.end())
    return i->second;

  // Subtract two for the initiating base pair.
  int length = en - st - 1;
  if (length < 3) return MAX_E;  // Disallowed by T04.
  energy_t energy = HairpinInitiation(length);

  // T04 says hairpin loops with all C bases inside them are treated specially.
  bool all_c = true;
  for (int i = st + 1; i <= en - 1; ++i) {
    if (r[i] != C) all_c = false;
  }

  if (length == 3) return energy + hairpin_c3_loop;
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
energy_t BulgeEnergy(int outer_st, int outer_en, int inner_st, int inner_en) {
  assert(
      inner_st > outer_st && inner_en < outer_en &&
      (outer_en - inner_en == 1 || inner_st - outer_st == 1) &&
      (outer_en - inner_en >= 2 || inner_st - outer_st >= 2));
  int length = std::max(inner_st - outer_st, outer_en - inner_en) - 1;
  energy_t energy = BulgeInitiation(length);
  // N.B. AU/GU penalty applied in ComputeEnergyInternal for bulges of length > 1.
  if (length > 1) return energy;
  // Stacking energy.
  energy += stacking_e[r[outer_st]][r[inner_st]][r[inner_en]][r[inner_st]];
  int unpaired = outer_st + 1;
  if (outer_st + 1 == inner_st) unpaired = inner_en + 1;
  // Special C bulge.
  if (r[unpaired] == C && (r[unpaired - 1] == C || r[unpaired] + 1 == C))
    energy += bulge_special_c;
#if BULGE_LOOP_STATES
  // TODO. This is not implemented and potentially experimentally not super solid.
#endif
  return energy;
}

energy_t ComputeEnergyInternal(int st, int en) {
  assert(en >= st - 1);
  // If () or .. or . or empty string, 0 energy.
  if (en - st <= 1) return 0;
  energy_t energy = 0;

  // Watson-Crick helix.
  if (st != -1) {
    if (p[st + 1] == en - 1) {
      printf("Adding Watson-Crick stacking energy %d %d %d %d: %d\n", st, st + 1, en - 1, en,
             stacking_e[r[st]][r[st + 1]][r[en - 1]][r[en]]);
      energy += stacking_e[r[st]][r[st + 1]][r[en - 1]][r[en]];
    } else if (IsUnorderedOf(r[st], r[en], G_b | A_b, U_b)) {
      // N.B. Double counting penalties for stacks of size 1 is intentional.
      printf("Applying ending AU/GU penalty %d %d\n", st, en);
      energy += AUGU_PENALTY;
    }
  }

  // Look for loops inside.
  int numloops = 0;
  int first_st = -1, first_en = -1;
  for (int i = st + 1; i <= en - 1; ++i) {
    int pair = p[i];
    assert(pair < en);
    if (pair != -1) {
      printf("Loop at %d %d\n", i, pair);
      // If not a continuation of stacking and not a size 1 bulge loop, apply AU/GU penalty.
      if ((i - st + en - pair > 3) && IsUnorderedOf(r[i], r[pair], G_b | A_b, U_b)) {
        printf("Applying opening AU/GU penalty %d %d\n", st + 1, en - 1);
        energy += AUGU_PENALTY;
      }
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
  } else if (numloops == 1) {
    energy += BulgeEnergy(st, en, first_st, first_en);
  }
  printf("Found %d loops\n", numloops);
  return energy;
}

// Destroys folded_rna_t.
energy_t ComputeEnergy(folded_rna_t& frna) {
  r = std::move(frna.r);
  p = std::move(frna.p);
  return ComputeEnergyInternal(-1, (int) r.size());
}

}
}
