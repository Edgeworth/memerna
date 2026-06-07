// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASE_PFN_PFN_EXTERIOR_H_
#define BACKENDS_BASE_PFN_PFN_EXTERIOR_H_

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/boltz_model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
void PfnExterior(const Primary& r, const BoltzModel& bm, PfnState& state, const PseudofreeCfgT& pf);

template <class PseudofreeCfgT>
inline BoltzEnergy PairedWithPf(const PseudofreeCfgT& pf, const BoltzDpArray& dp, int st, int en) {
  BoltzEnergy res = dp[st][en][PT_P];
  if (st > en) res *= Boltz(pf.Paired(en, st));
  return res;
}

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_PFN_PFN_EXTERIOR_H_
