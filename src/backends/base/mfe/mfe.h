// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASE_MFE_MFE_H_
#define BACKENDS_BASE_MFE_MFE_H_

#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

// Basic MFE folding.
void MfeDebug(const Primary& r, const Model::Ptr& m, DpState& state);

// Basic MFE folding.
void MfeOpt(const Primary& r, const Model::Ptr& m, DpState& state);

// Sparse folding.
void MfeSparseOpt(const Primary& r, const Model::Ptr& m, DpState& state);

// Sparse folding with Lyngso's algorithm.
void MfeLyngsoSparseOpt(const Primary& r, const Model::Ptr& m, DpState& state);

Energy MfeExterior(const Primary& r, const Model::Ptr& m, DpState& state);

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_MFE_MFE_H_
