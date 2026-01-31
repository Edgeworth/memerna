// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_STACK_MFE_MFE_EXTERIOR_H_
#define BACKENDS_STACK_MFE_MFE_EXTERIOR_H_

#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/dp.h"
#include "model/primary.h"

namespace mrna::md::stack {

Energy MfeExterior(const Primary& r, const Model::Ptr& m, DpState& state);

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_MFE_MFE_EXTERIOR_H_
