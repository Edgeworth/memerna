// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_PFN_PFN_EXTERIOR_H_
#define BACKENDS_CANON_PFN_PFN_EXTERIOR_H_

#include "api/energy/energy_cfg.h"
#include "backends/canon/dp.h"
#include "backends/canon/energy/boltz_model.h"
#include "model/primary.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
void PfnExterior(const Primary& r, const BoltzModel& bm, PfnState& state, const PseudofreeCfgT& pf);

template <class PseudofreeCfgT>
inline BoltzEnergy PairedWithPf(const PseudofreeCfgT& pf, const BoltzDpArray& dp, int st, int en) {
  BoltzEnergy res = dp[st][en][PT_P];
  if (st > en) res *= Boltz(pf.Paired(en, st));
  return res;
}

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_PFN_PFN_EXTERIOR_H_
