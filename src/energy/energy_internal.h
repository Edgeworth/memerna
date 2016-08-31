#ifndef MEMERNA_ENERGY_INTERNAL_H
#define MEMERNA_ENERGY_INTERNAL_H

#include <deque>
#include "common.h"

namespace memerna {
namespace energy {
namespace internal {

typedef std::deque<std::pair<Ctd, energy_t>> branch_ctd_t;

energy_t ComputeOptimalCtd(const secondary_t& secondary,
    const std::deque<int>& branches, int outer_idx,
    bool use_first_lu, branch_ctd_t* branch_ctds);

// Takes the list representation of ctds in |branch_ctds| for |branches| branches and
// writes it in per-base representation to |base_ctds|.
void BranchCtdsToBaseCtds(const secondary_t& secondary, const std::deque<int>& branches,
    const branch_ctd_t& branch_ctds, std::vector<Ctd>* base_ctds);
// Reads the per-base ctd representation from |base_ctds| for |branches| branches and
// writes it in list representation to |branch_ctds|.
energy_t BaseCtdsToBranchCtds(const secondary_t& secondary, const std::deque<int>& branches,
    const std::vector<Ctd>& base_ctds, branch_ctd_t* branch_ctds);

}
}
}

#endif  //MEMERNA_ENERGY_INTERNAL_H
