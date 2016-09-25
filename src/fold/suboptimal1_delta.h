#ifndef MEMERNA_SUBOPTIMAL1_DELTA_H
#define MEMERNA_SUBOPTIMAL1_DELTA_H

#include "common.h"
#include "fold/fold_internal.h"
#include "fold/suboptimal1_base.h"

namespace memerna {
namespace fold {
namespace internal {

class Suboptimal1Delta : public Suboptimal1Base {
public:
  Suboptimal1Delta(energy_t delta_) : Suboptimal1Base(delta_) {
    verify_expr(delta_ >= 0, "energy delta must be non-negative");
  }

  void Run(std::function<void(const computed_t&)> fn);

private:
  struct dfs_state_t {
    int idx;
    index_t expand;
    // Stores whether this node's |expand| was from |unexpanded| and needs to be replaced.
    bool should_unexpand;
  };

  std::stack<dfs_state_t> q;
  std::vector<index_t> unexpanded;
};
}
}
}

#endif  // MEMERNA_SUBOPTIMAL1_DELTA_H
