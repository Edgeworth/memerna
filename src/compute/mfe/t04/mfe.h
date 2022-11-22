#ifndef COMPUTE_MFE_T04_ALG_H_
#define COMPUTE_MFE_T04_ALG_H_

#include "compute/dp.h"
#include "compute/energy/energy.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::mfe::t04 {

// Basic MFE folding.
DpArray ComputeTablesSlowest(const Primary& r, const energy::EnergyModelPtr& em);

// Basic MFE folding.
DpArray ComputeTablesSlow(const Primary& r, const energy::EnergyModelPtr& em);

// Sparse folding.
DpArray ComputeTablesFastest(const Primary& r, const energy::EnergyModelPtr& em);

// Sparse folding with Lyngso's algorithm.
DpArray ComputeTablesLyngo(const Primary& r, const energy::EnergyModelPtr& em);

ExtArray ComputeExterior(const Primary& r, const energy::EnergyModel& em, const DpArray& dp);

}  // namespace mrna::mfe::t04

#endif  // COMPUTE_MFE_T04_ALG_H_
