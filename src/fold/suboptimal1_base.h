#ifndef MEMERNA_SUBOPTIMAL1_BASE_H
#define MEMERNA_SUBOPTIMAL1_BASE_H

#include <algorithm>
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

class Suboptimal1Base {
  Suboptimal1Base() = delete;

protected:
  Suboptimal1Base(energy_t delta_) : delta(delta_) {}

  const energy_t delta;

  const std::vector<expand_t>& GetExpansions(const index_t& to_expand) {
    auto iter = cache.find(to_expand);
    if (iter == cache.end()) {
      auto exps = GenerateExpansions(to_expand, delta);
      std::sort(exps.begin(), exps.end());
      auto res = cache.emplace(to_expand, std::move(exps));
      assert(res.second && res.first != cache.end());
      iter = res.first;
    }
    return iter->second;
  }

private:
  // This node is where we build intermediate results to be pushed onto the queue.
  std::unordered_map<index_t, std::vector<expand_t>> cache;
};

}
}
}

#endif  // MEMERNA_SUBOPTIMAL1_BASE_H
