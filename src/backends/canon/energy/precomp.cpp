// Copyright 2026 Eliot Courtney.
#include "backends/canon/energy/precomp.h"

#include "backends/common/static.h"

namespace mrna::md::base {

#define MRNA_INSTANTIATE(CFG, PF) template class Precomp<CFG, PF, mrna::md::canon::Model>;
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::base
