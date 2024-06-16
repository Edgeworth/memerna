// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_MFE_MFE_H_
#define BACKENDS_BASEOPT_MFE_MFE_H_

#include "backends/baseopt/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

// Debug folding.
void MfeDebug(const Primary& r, const Model::Ptr& m, DpState& state);

// Basic MFE folding.
void MfeOpt(const Primary& r, const Model::Ptr& m, DpState& state);

// Sparse folding.
void MfeSparseOpt(const Primary& r, const Model::Ptr& m, DpState& state);

// Sparse folding with Lyngso's algorithm.
void MfeLyngsoSparseOpt(const Primary& r, const Model::Ptr& m, DpState& state);

Energy MfeExterior(const Primary& r, const Model::Ptr& m, DpState& state);

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_MFE_MFE_H_
