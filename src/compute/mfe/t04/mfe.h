// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_MFE_T04_MFE_H_
#define COMPUTE_MFE_T04_MFE_H_

#include "compute/dp.h"
#include "compute/energy/t04/model.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::mfe::t04 {

// Basic MFE folding.
DpArray MfeSlowest(const Primary& r, const erg::t04::ModelPtr& em);

// Basic MFE folding.
DpArray MfeSlow(const Primary& r, const erg::t04::ModelPtr& em);

// Sparse folding.
DpArray MfeFastest(const Primary& r, const erg::t04::ModelPtr& em);

// Sparse folding with Lyngso's algorithm.
DpArray MfeLyngso(const Primary& r, const erg::t04::ModelPtr& em);

ExtArray MfeExterior(const Primary& r, const erg::t04::ModelPtr& em, const DpArray& dp);

}  // namespace mrna::mfe::t04

#endif  // COMPUTE_MFE_T04_MFE_H_
