// Copyright 2023 Eliot Courtney.
#ifndef COMPUTE_MFE_T22_DP_H_
#define COMPUTE_MFE_T22_DP_H_

#include <cassert>

#include "compute/mfe/t04/dp.h"
#include "model/base.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "util/array.h"

namespace mrna::mfe::t22 {

struct DpState {
  // T04 state is reused.
  t04::DpState dp;
  // DP for stack length.
  Array3D stack;
};

}  // namespace mrna::mfe::t22

#endif  // COMPUTE_MFE_T22_DP_H_
