// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_MFE_GLOBALS_H_
#define COMPUTE_MFE_GLOBALS_H_

#include <string>
#include <vector>

#include "common.h"
#include "compute/constants.h"
#include "compute/dp.h"
#include "compute/energy/model.h"
#include "util/array.h"

namespace mrna::mfe {

namespace internal {

extern std::vector<int> gp;
extern std::vector<Ctd> gctd;
extern std::string grep;
extern energy_t genergy;
extern array3d_t<energy_t, DP_SIZE> gdp;
extern array2d_t<energy_t, EXT_SIZE> gext;

}  // namespace internal

void SetMfeGlobalState(const primary_t& r, const energy::EnergyModel& em);

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_GLOBALS_H_
