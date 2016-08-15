#ifndef MEMERNA_TRACEBACK_H
#define MEMERNA_TRACEBACK_H

#include <stack>
#include "common.h"
#include "array.h"
#include "fold/fold.h"

namespace memerna {
namespace fold {

energy_t TraceExterior(const array3d_t<energy_t, DP_SIZE>& arr, std::stack<std::tuple<int, int, int>>& q);

void TraceStructure(const array3d_t<energy_t, DP_SIZE>& arr, std::stack<std::tuple<int, int, int>>& q);

}
}

#endif //MEMERNA_TRACEBACK_H
