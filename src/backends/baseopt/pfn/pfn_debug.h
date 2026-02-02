// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_PFN_PFN_DEBUG_H_
#define BACKENDS_BASEOPT_PFN_PFN_DEBUG_H_

#include <string>

#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "backends/baseopt/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/pfn.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

class PfnDebug {
 public:
  static bool IsSupported(
      const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

  static PfnTables Run(const Primary& r, const Model::Ptr& m, erg::EnergyCfg cfg, PfnState& state,
      const erg::PseudofreeCfg& pf);
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_PFN_PFN_DEBUG_H_
