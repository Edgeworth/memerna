// Copyright 2022 Eliot Courtney.
#include "backends/base/energy/precomp.h"

#include "backends/common/static.h"

namespace mrna::md::base {

#define MRNA_INSTANTIATE(CFG, PF) template class Precomp<CFG, PF>;
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::base
