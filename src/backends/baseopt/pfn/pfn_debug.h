// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_PFN_PFN_DEBUG_H_
#define BACKENDS_BASEOPT_PFN_PFN_DEBUG_H_

#include "backends/baseopt/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/pfn.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

class PfnDebug {
 public:
  static PfnTables Run(const Primary& r, const Model::Ptr& initial_m, PfnState& state);
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_PFN_PFN_DEBUG_H_
