#ifndef MEMERNA_CTDS_H
#define MEMERNA_CTDS_H

#include <deque>
#include "common.h"

namespace memerna {
namespace energy {

typedef std::deque<std::pair<Ctd, energy_t>> branch_ctd_t;

energy_t ComputeOptimalCtd(const std::deque<int>& branches, int outer_idx,
    bool use_first_lu, branch_ctd_t* ctd_energies);

void WriteCtdsForBranches(const std::deque<int>& branches, const branch_ctd_t& ctd_energies);
energy_t ReadCtdsForBranches(const std::deque<int>& branches, branch_ctd_t* ctd_energies);

}
}

#endif  //MEMERNA_CTDS_H
