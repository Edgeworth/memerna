// Copyright 2022 Eliot Courtney.
#ifndef MODELS_T04_MFE_MFE_H_
#define MODELS_T04_MFE_MFE_H_

#include "model/constants.h"
#include "model/primary.h"
#include "models/t04/energy/model.h"
#include "models/t04/mfe/dp.h"

namespace mrna::md::t04 {

// Basic MFE folding.
void MfeDebug(const Primary& r, const Model::Ptr& em, DpState& state);

// Basic MFE folding.
void MfeOpt(const Primary& r, const Model::Ptr& em, DpState& state);

// Sparse folding.
void MfeSparseOpt(const Primary& r, const Model::Ptr& em, DpState& state);

// Sparse folding with Lyngso's algorithm.
void MfeLyngsoSparseOpt(const Primary& r, const Model::Ptr& em, DpState& state);

Energy MfeExterior(const Primary& r, const Model::Ptr& em, DpState& state);

}  // namespace mrna::md::t04

#endif  // MODELS_T04_MFE_MFE_H_
