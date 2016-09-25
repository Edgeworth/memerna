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
  Suboptimal1Delta(energy_t delta_, int num) :
      Suboptimal1Base(delta_ == -1 ? CAP_E : delta_),
      max_structures(num == -1 ? MAX_STRUCTURES : num) {}

  int Run(std::function<void(const computed_t&)> fn, bool sorted);

private:
  struct dfs_state_t {
    int idx;
    index_t expand;
    // Stores whether this node's |expand| was from |unexpanded| and needs to be replaced.
    bool should_unexpand;
  };

  constexpr static int MAX_STRUCTURES = std::numeric_limits<int>::max() / 4;

  const int max_structures;
  std::stack<dfs_state_t> q;
  std::vector<index_t> unexpanded;

  std::pair<int, int> RunInternal(std::function<void(const computed_t&)> fn,
    energy_t cur_delta, bool exact_energy, int structure_limit);
};
}
}
}

#endif  // MEMERNA_SUBOPTIMAL1_DELTA_H
