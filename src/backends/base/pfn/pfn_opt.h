// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASE_PFN_PFN_OPT_H_
#define BACKENDS_BASE_PFN_PFN_OPT_H_

#include <string>

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/boltz_model.h"
#include "backends/common/base/dp.h"
#include "model/pfn.h"
#include "model/primary.h"

namespace mrna::md::base {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
class PfnOpt {
 public:
  static bool IsSupported(const PseudofreeCfgT& pf, std::string* reason = nullptr);

  static PfnTables Run(
      const Primary& r, const BoltzModel::Ptr& bm, PfnState& state, const PseudofreeCfgT& pf);
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_PFN_PFN_OPT_H_
