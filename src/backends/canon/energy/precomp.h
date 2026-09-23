// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_ENERGY_PRECOMP_H_
#define BACKENDS_CANON_ENERGY_PRECOMP_H_

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/precomp.h"
#include "backends/canon/energy/model.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
using Precomp = base::Precomp<CFG, PseudofreeCfgT, Model>;

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_ENERGY_PRECOMP_H_
