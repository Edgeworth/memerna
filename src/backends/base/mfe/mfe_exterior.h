// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_BASE_MFE_MFE_EXTERIOR_H_
#define BACKENDS_BASE_MFE_MFE_EXTERIOR_H_

#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

Energy MfeExterior(const Primary& r, const Model::Ptr& m, DpState& state);

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_MFE_MFE_EXTERIOR_H_
