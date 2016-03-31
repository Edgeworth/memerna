#include <cstdio>
#include <cmath>
#include "energy.h"
#include "globals.h"

namespace memerna {
namespace energy {
namespace internal {

energy_t HairpinInitiation(int n) {
  assert(n > 0);
  if (n < INITIATION_CACHE_SZ) return hairpin_init[n - 1];
  static_assert(INITIATION_CACHE_SZ > 9, "Need initiation values for up to 9.");
  // Formula: G_init(9) + 1.75 * R * T * ln(n / 9).
  return static_cast<energy_t>(round(hairpin_init[8] + 1.75 * R * T * log(n / 9.0)));
}

// Indices are inclusive, include the initiating base pair.
// Rules for hairpin energy:
// 1. Check a lookup table for a hardcoded value.
// 2. Loops with less than 3 bases on the inside are disallowed by T04.
// 3. If it has 3 bases on the inside, return G_init + penalty for all C bases.
// 3.1 The penalty for all C bases for size 3 hairpin loops is special cased.
// 4. If it has more than 3 bases on the inside... TODO
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
  base_t left = r[st+1], right = r[en-1];
  energy += terminal_e[r[st]][left][right][r[en]];
  if (IsPairOf(left, right, U_b, U_b) || IsPairOf(left, right, G_b, A_b))
    energy += hairpin_uu_ga_first_mismatch;
  if (IsPairOf(left, right, G_b, G_b))
    energy += hairpin_gg_first_mismatch;
  if (all_c)
    energy += hairpin_all_c_a * length + hairpin_all_c_b;

  assert(st >= 2);
  if (IsPairOf(r[st], r[en], G_b, U_b) && r[st-1] == G && r[st-2] == G)
    energy += hairpin_special_gu_closure;
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
  for (int i = st + 1; i <= en - 1; ++i) {
    int pair = p[i];
    assert(pair < en);
    if (pair != -1) {
      printf("Loop at %d %d\n", i, pair);
      // If not a continuation of stacking, apply AU/GU penalty.
      if ((i != st + 1 || pair != en - 1) && IsUnorderedOf(r[i], r[pair], G_b | A_b, U_b)) {
        printf("Applying opening AU/GU penalty %d %d\n", st + 1, en - 1);
        energy += AUGU_PENALTY;
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
  }
  printf("Found %d loops\n", numloops);
  return energy;
}
}

// Destroys folded_rna_t.
energy_t ComputeEnergy(folded_rna_t& frna) {
  r = std::move(frna.r);
  p = std::move(frna.p);
  return internal::ComputeEnergyInternal(-1, (int) r.size());
}

}
}
