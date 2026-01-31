// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_MFE_MFE_OPT_H_
#define BACKENDS_BASEOPT_MFE_MFE_OPT_H_

#include "backends/baseopt/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

// Basic MFE folding.
class MfeOpt {
 public:
  static void Run(const Primary& r, const Model::Ptr& m, DpState& state);
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_MFE_MFE_OPT_H_
