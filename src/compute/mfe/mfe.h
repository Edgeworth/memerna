// Copyright 2016 E.
#ifndef COMPUTE_MFE_MFE_H_
#define COMPUTE_MFE_MFE_H_

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::mfe {

struct MfeResult {
  DpArray dp;
  ExtArray ext;
  Energy energy = ZERO_E;
};

// TODO(3): Move this?
struct Cand {
  Energy energy;
  int idx;
};

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_MFE_H_
