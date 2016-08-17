#ifndef MEMERNA_ENERGY_H
#define MEMERNA_ENERGY_H

#include <utility>
#include <string>
#include <deque>
#include <memory>
#include "common.h"
#include "base.h"
#include "globals.h"

namespace memerna {

namespace structure {
class Structure;
}

namespace energy {

inline energy_t AuGuPenalty(int st, int en) {
  return IsAuGu(r[st], r[en]) ? g_augu_penalty : 0;
}

energy_t HairpinInitiation(int n);

energy_t Hairpin(int st, int en, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t BulgeInitiation(int n);

energy_t Bulge(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t InternalLoopInitiation(int n);

energy_t InternalLoop(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t TwoLoop(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t MultiloopInitiation(int num_branches);

// We use the normal terminal mismatch parameters for the mismatch that is on the continuous part of the
// RNA. The stacking for the non-continuous part is set to be an arbitrary given number. There are two possible
// orientations, since the base involved in the terminal mismatch could come from either side.
// ... _ _ _ _ ...
// ...|_|  _|_|...
//      | |
// Rules for mismatch mediated coaxial stacking:
//    1. A terminal mismatch is formed around the branch being straddled.
//    2. An arbitrary bonus is added.
//    2. An arbitrary bonus is added if the mismatch is Watson-Crick or GU.
inline energy_t MismatchCoaxial(base_t fiveTop, base_t mismatch_top, base_t mismatch_bot, base_t threeBot) {
  energy_t coax = g_terminal[fiveTop][mismatch_top][mismatch_bot][threeBot] + g_coax_mismatch_non_contiguous;
  if (IsWatsonCrick(mismatch_top, mismatch_bot))
    coax += g_coax_mismatch_wc_bonus;
  else if (IsGu(mismatch_top, mismatch_bot))
    coax += g_coax_mismatch_gu_bonus;
  return coax;
}


// Requires global variables r and p to be set. st and en are inclusive.
energy_t ComputeEnergy(std::unique_ptr<structure::Structure>* s = nullptr);

inline energy_t ComputeEnergy(const folded_rna_t& frna, std::unique_ptr<structure::Structure>* s = nullptr) {
  SetFoldedRna(frna);
  return ComputeEnergy(s);
}

}
}

#endif //MEMERNA_ENERGY_H
