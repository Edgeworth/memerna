#include <cstdio>
#include "energy.h"
#include "globals.h"

namespace memerna {
namespace energy {
namespace internal {

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

  // Bulge loop.
  if (numloops == 0) {

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
