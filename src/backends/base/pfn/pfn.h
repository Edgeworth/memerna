// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASE_PFN_PFN_H_
#define BACKENDS_BASE_PFN_PFN_H_

#include "backends/base/energy/boltz_model.h"
#include "backends/common/base/dp.h"
#include "model/pfn.h"
#include "model/primary.h"

namespace mrna::md::base {

PfnTables PfnDebug(const Primary& r, const Model::Ptr& initial_m, PfnState& state);

PfnTables PfnOpt(const Primary& r, const BoltzModel::Ptr& bm, PfnState& state);

void PfnExterior(const Primary& r, const ModelBase& m, PfnState& state);

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_PFN_PFN_H_
