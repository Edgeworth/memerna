#ifndef MEMERNA_SUBOPTIMAL1_H
#define MEMERNA_SUBOPTIMAL1_H

#include "common.h"
#include "fold/fold_internal.h"

namespace memerna {
namespace fold {
namespace internal {

struct expand_t {
  expand_t() = delete;
  expand_t(energy_t energy_) : energy(energy_) {}
  expand_t(energy_t energy_, const index_t& to_expand_) : energy(energy_), to_expand(to_expand_) {}
  expand_t(energy_t energy_, const index_t& to_expand_, const index_t& unexpanded_)
      : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_) {}
  expand_t(energy_t energy_, const index_t& to_expand_, const ctd_idx_t& ctd0_)
      : energy(energy_), to_expand(to_expand_), ctd0(ctd0_) {}
  expand_t(energy_t energy_, const index_t& to_expand_, const index_t& unexpanded_,
      const ctd_idx_t& ctd0_)
      : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_), ctd0(ctd0_) {}
  expand_t(energy_t energy_, const index_t& to_expand_, const index_t& unexpanded_,
      const ctd_idx_t& ctd0_, const ctd_idx_t& ctd1_)
      : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_), ctd0(ctd0_), ctd1(ctd1_) {}
  energy_t energy;
  index_t to_expand, unexpanded;  // st is -1 if this does not exist
  ctd_idx_t ctd0, ctd1;

  // TODO return ones with less work first.
  bool operator<(const expand_t& o) const { return energy < o.energy; }
};

std::vector<expand_t> GenerateExpansions(const index_t& to_expand, energy_t delta);


class Suboptimal1 {
public:
  Suboptimal1(energy_t delta_, int num) :
      delta(delta_ == -1 ? CAP_E : delta_),
      max_structures(num == -1 ? MAX_STRUCTURES : num) {}

  int Run(std::function<void(const computed_t&)> fn, bool sorted);

private:
  constexpr static int MAX_STRUCTURES = std::numeric_limits<int>::max() / 4;

  struct dfs_state_t {
    int idx;
    index_t expand;
    // Stores whether this node's |expand| was from |unexpanded| and needs to be replaced.
    bool should_unexpand;
  };

  const energy_t delta;
  const int max_structures;
  // This node is where we build intermediate results to be pushed onto the queue.
  std::unordered_map<index_t, std::vector<expand_t>> cache;
  std::stack<dfs_state_t> q;
  std::vector<index_t> unexpanded;

  std::pair<int, int> RunInternal(std::function<void(const computed_t&)> fn,
    energy_t cur_delta, bool exact_energy, int structure_limit);

  const std::vector<expand_t>& GetExpansions(const index_t& to_expand);
};
}
}
}

#endif  // MEMERNA_SUBOPTIMAL1_H
