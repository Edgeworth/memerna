// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_STACK_MFE_MFE_OPT_H_
#define BACKENDS_STACK_MFE_MFE_OPT_H_

#include <string>

#include "api/energy/pseudofree_cfg.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/dp.h"
#include "model/primary.h"

namespace mrna::md::stack {

class MfeOpt {
 public:
  static bool IsSupported(
      const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

  static void Run(const Primary& r, const Model::Ptr& m, DpState& state, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf);
};

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_MFE_MFE_OPT_H_
