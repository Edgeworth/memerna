#ifndef MEMERNA_SUBOPTIMAL_H
#define MEMERNA_SUBOPTIMAL_H

#include <stack>
#include "common.h"
#include "array.h"
#include "fold/fold.h"

namespace memerna {
namespace fold {

std::vector<folded_rna_t> SuboptimalTraceback0(
    energy_t max_energy, int max_structures,
    const array3d_t<energy_t, DP_SIZE>& arr,
    const array2d_t<energy_t, EXT_SIZE>& exterior);

}
}

#endif   //MEMERNA_SUBOPTIMAL_H
