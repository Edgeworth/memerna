// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_STACK_MFE_MFE_DEBUG_H_
#define BACKENDS_STACK_MFE_MFE_DEBUG_H_

#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/dp.h"
#include "model/primary.h"

namespace mrna::md::stack {

class MfeDebug {
 public:
  static void Run(const Primary& r, const Model::Ptr& m, DpState& state);
};

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_MFE_MFE_DEBUG_H_
