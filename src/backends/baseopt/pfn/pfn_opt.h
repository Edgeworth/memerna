// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_PFN_PFN_OPT_H_
#define BACKENDS_BASEOPT_PFN_PFN_OPT_H_

#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "backends/baseopt/energy/boltz_model.h"
#include "backends/common/base/dp.h"
#include "model/pfn.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

class PfnOpt {
 public:
  static PfnTables Run(const Primary& r, const BoltzModel::Ptr& bm, erg::EnergyCfg cfg,
      PfnState& state, const erg::PseudofreeCfg& pf);
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_PFN_PFN_OPT_H_
