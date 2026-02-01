// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASE_MFE_MFE_LYNGSO_SPARSE_OPT_H_
#define BACKENDS_BASE_MFE_MFE_LYNGSO_SPARSE_OPT_H_

#include <string>

#include "api/energy/pseudofree_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

// Sparse folding with Lyngso's algorithm.
class MfeLyngsoSparseOpt {
 public:
  static bool IsSupported(
      const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

  static void Run(const Primary& r, const Model::Ptr& m, DpState& state, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf);
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_MFE_MFE_LYNGSO_SPARSE_OPT_H_
