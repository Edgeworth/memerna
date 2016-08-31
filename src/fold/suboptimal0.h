#ifndef MEMERNA_SUBOPTIMAL0_H
#define MEMERNA_SUBOPTIMAL0_H

#include <stack>
#include "common.h"
#include "array.h"
#include "fold/fold.h"

namespace memerna {
namespace fold {

std::vector<computed_t> SuboptimalTraceback0(
    energy_t max_energy, int max_structures,
    const primary_t& r,
    const array3d_t<energy_t, DP_SIZE>& arr,
    const array2d_t<energy_t, EXT_SIZE>& exterior);

}
}

#endif   //MEMERNA_SUBOPTIMAL0_H
