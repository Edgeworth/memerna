// Copyright 2023 Eliot Courtney.

#ifndef COMPUTE_MFE_T22_MFE_H_
#define COMPUTE_MFE_T22_MFE_H_

#include "compute/dp.h"
#include "compute/energy/t22/model.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::mfe::t22 {

DpArray MfeSlowest(const Primary& r, const erg::t22::ModelPtr& em);

ExtArray MfeExterior(const Primary& r, const erg::t22::ModelPtr& em, const DpArray& dp);

}  // namespace mrna::mfe::t22

#endif  // COMPUTE_MFE_T22_MFE_H_
