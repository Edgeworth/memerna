#ifndef MEMERNA_SLOW_FOLD_H
#define MEMERNA_SLOW_FOLD_H

#include "common.h"
#include "array.h"

namespace memerna {
namespace fold {

array3d_t<energy_t, DP_SIZE> ComputeTables_Slow();

}
}

#endif //MEMERNA_SLOW_FOLD_H
