// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_ENERGY_BOLTZ_PRECOMP_H_
#define BACKENDS_CANON_ENERGY_BOLTZ_PRECOMP_H_

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/boltz_precomp.h"
#include "backends/canon/energy/boltz_model.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
using BoltzPrecomp = base::BoltzPrecomp<CFG, PseudofreeCfgT, BoltzModel>;

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_ENERGY_BOLTZ_PRECOMP_H_
