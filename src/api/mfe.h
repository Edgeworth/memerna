// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_MFE_MFE_H_
#define COMPUTE_MFE_MFE_H_

#include "api/energy/model.h"
#include "model/constants.h"
#include "model/primary.h"
#include "models/t04/mfe/dp.h"

namespace mrna::mfe {

using DpState = std::variant<md::t04::DpState>;

struct MfeResult {
  DpState dp;
  Energy energy = ZERO_E;
};

}  // namespace mrna::mfe

#endif  // COMPUTE_MFE_MFE_H_
