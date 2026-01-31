// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASE_PFN_PFN_EXTERIOR_H_
#define BACKENDS_BASE_PFN_PFN_EXTERIOR_H_

#include "backends/base/energy/boltz_model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

void PfnExterior(const Primary& r, const Model& m, PfnState& state);

inline BoltzEnergy PairedWithPf(const Model::Ptr& m, const BoltzDpArray& dp, int st, int en) {
  BoltzEnergy res = dp[st][en][PT_P];
  if (st > en) res *= m->pf.Paired(en, st).Boltz();
  return res;
}

inline BoltzEnergy PairedWithPf(const BoltzModel::Ptr& bm, const BoltzDpArray& dp, int st, int en) {
  BoltzEnergy res = dp[st][en][PT_P];
  if (st > en) res *= bm->pf.Paired(en, st);
  return res;
}

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_PFN_PFN_EXTERIOR_H_
