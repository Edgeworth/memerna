// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_MFE_MFE_EXTERIOR_H_
#define BACKENDS_CANON_MFE_MFE_EXTERIOR_H_

#include "backends/canon/dp.h"
#include "backends/canon/energy/model.h"
#include "model/primary.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
Energy MfeExterior(const Primary& r, const Model::Ptr& m, DpState& state, const PseudofreeCfgT& pf);

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_MFE_MFE_EXTERIOR_H_
