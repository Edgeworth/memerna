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
  return IsAuGu(r[st], r[en]) ? augu_penalty : 0;
}

energy_t HairpinInitiation(int n);

energy_t HairpinEnergy(int st, int en, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t BulgeInitiation(int n);

energy_t BulgeEnergy(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t InternalLoopInitiation(int n);

energy_t InternalLoopEnergy(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t TwoLoop(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t MultiloopInitiation(int num_branches);

energy_t MismatchMediatedCoaxialEnergy(
    base_t fiveTop, base_t mismatch_top, base_t mismatch_bot, base_t threeBot);

// Requires global variables r and p to be set. st and en are inclusive.
energy_t ComputeEnergy(std::unique_ptr<structure::Structure>* s = nullptr);

inline energy_t ComputeEnergy(const folded_rna_t& frna, std::unique_ptr<structure::Structure>* s = nullptr) {
  SetFoldedRna(frna);
  return ComputeEnergy(s);
}

}
}

#endif //MEMERNA_ENERGY_H
