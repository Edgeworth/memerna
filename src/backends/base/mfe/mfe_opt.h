// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASE_MFE_MFE_OPT_H_
#define BACKENDS_BASE_MFE_MFE_OPT_H_

#include <string>

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

// Basic MFE folding.
template <erg::EnergyCfg CFG, class PseudofreeCfgT>
class MfeOpt {
 public:
  static bool IsSupported(const PseudofreeCfgT& pf, std::string* reason = nullptr);

  static void Run(const Primary& r, const Model::Ptr& m, DpState& state, const PseudofreeCfgT& pf);
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_MFE_MFE_OPT_H_
