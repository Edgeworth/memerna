#ifndef MEMERNA_BRUTE_FOLD_H
#define MEMERNA_BRUTE_FOLD_H

#include "common.h"
#include "energy/energy_model.h"

namespace memerna {
namespace fold {

computed_t FoldBruteForce(const primary_t& r, const energy::EnergyModel& em);

}
}
#endif  // MEMERNA_BRUTE_FOLD_H
