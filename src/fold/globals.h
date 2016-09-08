#ifndef MEMERNA_FOLD_GLOBALS_H
#define MEMERNA_FOLD_GLOBALS_H

#include "common.h"
#include "array.h"
#include "energy/energy_model.h"
#include "fold/fold_constants.h"
#include "fold/precomp.h"

namespace memerna {
namespace fold {
namespace internal {

extern primary_t gr;
extern std::vector<int> gp;
extern std::vector<Ctd> gctd;
extern energy_t genergy;
extern energy::EnergyModel gem;
extern precomp_t gpc;
extern array3d_t<energy_t, DP_SIZE> gdp;
extern array2d_t<energy_t, EXT_SIZE> gext;

void SetGlobalState(const primary_t& r, const energy::EnergyModel& em);

}
}
}

#endif  // MEMERNA_FOLD_GLOBALS_H
