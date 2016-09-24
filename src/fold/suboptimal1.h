#ifndef MEMERNA_SUBOPTIMAL1_H
#define MEMERNA_SUBOPTIMAL1_H

#include <algorithm>
#include <set>
#include "common.h"
#include "energy/structure.h"
#include "fold/fold_internal.h"
#include "fold/suboptimal1_delta.h"
#include "fold/suboptimal1_num.h"
#include "parsing.h"

namespace memerna {
namespace fold {
namespace internal {

class Suboptimal1 {
public:
  Suboptimal1(energy_t delta_, int num_) : delta(delta_), num(num_) {}

  void Run(std::function<void(const computed_t&)> fn) {
    verify_expr(
        gr.size() < std::numeric_limits<int16_t>::max(), "RNA too long for suboptimal folding");
    if (num == -1) Suboptimal1Delta(delta).Run(fn);
    else Suboptimal1Num(delta, num).Run(fn);
  }

private:
  const energy_t delta;
  const int num;
};
}
}
}

#endif  // MEMERNA_SUBOPTIMAL1_H
