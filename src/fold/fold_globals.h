// Copyright 2016 E.
#ifndef FOLD_FOLD_GLOBALS_H_
#define FOLD_FOLD_GLOBALS_H_

#include <string>
#include <vector>

#include "array.h"
#include "common.h"
#include "energy/energy_model.h"
#include "fold/fold_constants.h"

namespace mrna {
namespace fold {
namespace internal {

extern std::vector<int> gp;
extern std::vector<Ctd> gctd;
extern std::string grep;
extern energy_t genergy;
extern array3d_t<energy_t, DP_SIZE> gdp;
extern array2d_t<energy_t, EXT_SIZE> gext;

}  // namespace internal

void SetFoldGlobalState(const primary_t& r, const energy::EnergyModel& em);

}  // namespace fold
}  // namespace mrna

#endif  // FOLD_FOLD_GLOBALS_H_
