// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_MFE_MFE_H_
#define COMPUTE_MFE_MFE_H_

#include "compute/energy/model.h"
#include "compute/mfe/t04/dp.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::mfe {

using DpState = std::variant<t04::DpState>;

struct MfeResult {
  DpState dp;
  Energy energy = ZERO_E;
};

// TODO(3): Move this?
struct Cand {
  Energy energy;
  int idx;
};

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_MFE_H_
