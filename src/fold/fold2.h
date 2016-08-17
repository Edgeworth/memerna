#ifndef MEMERNA_FOLD2_H
#define MEMERNA_FOLD2_H

#include "common.h"
#include "fold.h"
#include "array.h"

namespace memerna {
namespace fold {

// Note that this doesn't completely follow the loaded energy model.
array3d_t<energy_t, DP_SIZE> ComputeTables2();

}
}

#endif //MEMERNA_FOLD2_H
