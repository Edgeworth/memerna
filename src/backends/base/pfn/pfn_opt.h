// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASE_PFN_PFN_OPT_H_
#define BACKENDS_BASE_PFN_PFN_OPT_H_

#include "api/energy/pseudofree_cfg.h"
#include "backends/base/energy/boltz_model.h"
#include "backends/common/base/dp.h"
#include "model/pfn.h"
#include "model/primary.h"

namespace mrna::md::base {

class PfnOpt {
 public:
  static PfnTables Run(
      const Primary& r, const BoltzModel::Ptr& bm, PfnState& state, const erg::PseudofreeCfg& pf);
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_PFN_PFN_OPT_H_
