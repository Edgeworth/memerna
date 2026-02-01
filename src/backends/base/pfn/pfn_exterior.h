// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASE_PFN_PFN_EXTERIOR_H_
#define BACKENDS_BASE_PFN_PFN_EXTERIOR_H_

#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

void PfnExterior(const Primary& r, const Model& m, erg::EnergyCfg cfg, PfnState& state,
    const erg::PseudofreeCfg& pf);

inline BoltzEnergy PairedWithPf(
    const erg::PseudofreeCfg& pf, const BoltzDpArray& dp, int st, int en) {
  BoltzEnergy res = dp[st][en][PT_P];
  if (st > en) res *= pf.Paired(en, st).Boltz();
  return res;
}

inline BoltzEnergy PairedWithPf(
    const erg::BoltzPseudofreeCfg& bpf, const BoltzDpArray& dp, int st, int en) {
  BoltzEnergy res = dp[st][en][PT_P];
  if (st > en) res *= bpf.Paired(en, st);
  return res;
}

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_PFN_PFN_EXTERIOR_H_
