#ifndef MEMERNA_ENERGY_H
#define MEMERNA_ENERGY_H

#include <utility>
#include <string>
#include <deque>
#include <memory>
#include "common.h"
#include "base.h"

namespace memerna {

namespace structure {
class Structure;
}

namespace energy {

const energy_t MAX_E = 1000000000;
const energy_t AUGU_PENALTY = 5;
// N.B. This is for kcal/mol so it's not 8.315.
const double R = 1.985877534e-3;
// This is 37 degrees Celsius. Changing this is not a good idea.
const double T = 310.15;

energy_t HairpinInitiation(int n);

energy_t HairpinEnergy(int st, int en, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t BulgeInitiation(int n);

energy_t BulgeEnergy(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t InternalLoopInitiation(int n);

energy_t InternalLoopEnergy(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t TwoLoop(int ost, int oen, int ist, int ien, std::unique_ptr<structure::Structure>* s = nullptr);

energy_t MultiloopInitiation(int num_unpaired, int num_branches);

energy_t MultiloopT99Initiation(int num_unpaired, int num_branches);

energy_t MultiloopHackInitiation(int num_branches);

// This function is not re-entrant.
energy_t ComputeEnergy(const folded_rna_t& frna, std::unique_ptr<structure::Structure>* s = nullptr);

}
}

#endif //MEMERNA_ENERGY_H
