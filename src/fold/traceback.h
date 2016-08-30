#ifndef MEMERNA_TRACEBACK_H
#define MEMERNA_TRACEBACK_H

#include <stack>
#include "common.h"
#include "array.h"
#include "fold/fold.h"

namespace memerna {
namespace fold {

typedef std::stack<std::tuple<int, int, int>> traceback_stack_t;

array2d_t<energy_t, EXT_SIZE> TraceExterior(
    const array3d_t<energy_t, DP_SIZE>& arr, traceback_stack_t& q);
computed_t TraceStructure(const array3d_t<energy_t, DP_SIZE>& arr,
    const array2d_t<energy_t, EXT_SIZE>& exterior, traceback_stack_t& q);

}
}

#endif //MEMERNA_TRACEBACK_H
