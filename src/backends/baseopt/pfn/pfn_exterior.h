// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_PFN_PFN_EXTERIOR_H_
#define BACKENDS_BASEOPT_PFN_PFN_EXTERIOR_H_

#include "api/energy/energy_cfg.h"
#include "backends/baseopt/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

void PfnExterior(const Primary& r, const Model& m, erg::EnergyCfg cfg, PfnState& state);

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_PFN_PFN_EXTERIOR_H_
