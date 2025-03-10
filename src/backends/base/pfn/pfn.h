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

void PfnExterior(const Primary& r, const Model& m, PfnState& state);

inline BoltzEnergy PairedWithPf(const Model::Ptr& m, const BoltzDpArray& dp, int st, int en) {
  BoltzEnergy res = dp[st][en][PT_P];
  if (st > en) res *= m->pf.Paired(en, st).Boltz();
  return res;
}

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_PFN_PFN_H_
