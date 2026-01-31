// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_MFE_MFE_LYNGSO_SPARSE_OPT_H_
#define BACKENDS_BASEOPT_MFE_MFE_LYNGSO_SPARSE_OPT_H_

#include "api/energy/pseudofree_cfg.h"
#include "backends/baseopt/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

// Sparse folding with Lyngso's algorithm.
class MfeLyngsoSparseOpt {
 public:
  static void Run(
      const Primary& r, const Model::Ptr& m, DpState& state, const erg::PseudofreeCfg& pf);
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_MFE_MFE_LYNGSO_SPARSE_OPT_H_
