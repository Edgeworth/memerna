#ifndef MEMERNA_ENERGY_H
#define MEMERNA_ENERGY_H

#include <utility>
#include "common.h"
#include "base.h"

namespace memerna {
namespace energy {

typedef int32_t energy_t;
const energy_t MAX_E = 1000000000;
const energy_t AUGU_PENALTY = 5;
// N.B. This is for kcal/mol so it's not 8.315.
const double R = 1.985877534e-3;
// This is 37 degrees Celsius. Changing this is not a good idea.
const double T = 310.15;

energy_t HairpinInitiation(int n);

energy_t HairpinEnergy(int st, int en);

energy_t BulgeInitiation(int n);

energy_t BulgeEnergy(int outer_st, int outer_en, int inner_st, int inner_en);

energy_t ComputeEnergy(folded_rna_t& frna);
}
}

#endif //MEMERNA_ENERGY_H
