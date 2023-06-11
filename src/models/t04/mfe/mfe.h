// Copyright 2022 Eliot Courtney.
#ifndef MODELS_T04_MFE_MFE_H_
#define MODELS_T04_MFE_MFE_H_

#include "model/constants.h"
#include "model/primary.h"
#include "models/t04/energy/model.h"
#include "models/t04/mfe/dp.h"

namespace mrna::md::t04 {

// Basic MFE folding.
void MfeSlowest(const Primary& r, const Model::Ptr& em, DpState& state);

// Basic MFE folding.
void MfeSlow(const Primary& r, const Model::Ptr& em, DpState& state);

// Sparse folding.
void MfeFastest(const Primary& r, const Model::Ptr& em, DpState& state);

// Sparse folding with Lyngso's algorithm.
void MfeLyngso(const Primary& r, const Model::Ptr& em, DpState& state);

Energy MfeExterior(const Primary& r, const Model::Ptr& em, DpState& state);

}  // namespace mrna::md::t04

#endif  // MODELS_T04_MFE_MFE_H_
