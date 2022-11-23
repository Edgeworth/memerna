#ifndef COMPUTE_ENERGY_T04_BRANCH_H_
#define COMPUTE_ENERGY_T04_BRANCH_H_

#include <deque>
#include <utility>
#include <vector>

#include "compute/energy/branch.h"
#include "compute/energy/energy.h"
#include "compute/energy/t04/model.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::erg::t04 {

Energy ComputeOptimalCtds(const Model& em, const Primary& r, const Secondary& s,
    const std::deque<int>& branches, bool use_first_lu, BranchCtd* branch_ctd);

// Reads the per-base ctd representation from |ctd| for |branches| branches and
// writes it in branch representation to |branch_ctd|.
Energy AddBaseCtdsToBranchCtds(const Model& em, const Primary& r, const Secondary& s,
    const Ctds& ctd, const std::deque<int>& branches, BranchCtd* branch_ctd);

}  // namespace mrna::erg::t04

#endif  // COMPUTE_ENERGY_T04_BRANCH_H_
