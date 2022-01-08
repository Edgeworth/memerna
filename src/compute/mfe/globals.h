// Copyright 2016 E.
#ifndef COMPUTE_MFE_GLOBALS_H_
#define COMPUTE_MFE_GLOBALS_H_

#include <string>
#include <vector>

#include "compute/constants.h"
#include "compute/dp.h"
#include "compute/energy/model.h"
#include "util/array.h"

namespace mrna::mfe {

namespace internal {

extern std::vector<int> gp;
extern std::vector<Ctd> gctd;
extern std::string grep;
extern Energy genergy;
extern Array3D<Energy, DP_SIZE> gdp;
extern Array2D<Energy, EXT_SIZE> gext;

}  // namespace internal

void SetMfeGlobalState(const Primary& r);

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_GLOBALS_H_
