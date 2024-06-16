// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_PFN_PFN_H_
#define BACKENDS_BASEOPT_PFN_PFN_H_

#include "backends/baseopt/energy/boltz_model.h"
#include "backends/baseopt/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/pfn.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

PfnTables PfnDebug(const Primary& r, const Model::Ptr& initial_m, PfnState& state);

PfnTables PfnOpt(const Primary& r, const BoltzModel::Ptr& bm, PfnState& state);

void PfnExterior(const Primary& r, const Model& m, PfnState& state);

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_PFN_PFN_H_
