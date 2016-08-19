#ifndef MEMERNA_FOLD3_H
#define MEMERNA_FOLD3_H

#include "common.h"
#include "fold/fold.h"
#include "array.h"

namespace memerna {
namespace fold {

// Note that this doesn't completely follow the loaded energy model.
array3d_t<energy_t, DP_SIZE> ComputeTables3();

}
}

#endif //MEMERNA_FOLD3_H
