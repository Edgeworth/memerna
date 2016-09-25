#ifndef MEMERNA_SUBOPTIMAL1_NUM_H
#define MEMERNA_SUBOPTIMAL1_NUM_H

#include <set>
#include "common.h"
#include "fold/fold_internal.h"
#include "fold/suboptimal1_base.h"

namespace memerna {
namespace fold {
namespace internal {

class Suboptimal1Num : public Suboptimal1Base {
public:
  Suboptimal1Num(energy_t delta_, int num)
      : Suboptimal1Base(delta_ == -1 ? CAP_E : delta_),
        max_structures(num == -1 ? std::numeric_limits<int>::max() / 4 : num),
        q([this](int a, int b) { return NodeComparator(a, b); }) {
    verify_expr(max_structures > 0, "must request at least one structure");
  }

  void Run(std::function<void(const computed_t&)> fn);

private:
  struct num_node_t {
    num_node_t() = delete;
    num_node_t(const expand_t& exp_)
        : exp(exp_), child_count(0), parent(-1), exp_st(0), exp_en(-1), cur_anc(0) {}
    expand_t exp;
    // See Suboptimal1Delta for explanation.
    int child_count, parent, exp_st, exp_en, cur_anc;
  };

  const int max_structures;
  std::vector<num_node_t> nodes;
  std::vector<int> finished;  // TODO get rid of this, don't need it with new reconstruction method.
  std::vector<int> free_space;                     // Free space from stuff that has been GC'd.
  std::set<int, std::function<bool(int, int)>> q;  // Queue of nodes ready to be expanded.

  bool NodeComparator(int a, int b) const {
    if (nodes[a].exp.energy != nodes[b].exp.energy)
      return nodes[a].exp.energy < nodes[b].exp.energy;
    return a < b;
  }

  computed_t ReconstructComputed(int node_idx);
  std::vector<computed_t> ConstructComputedsFromNodes();
  bool InsertQ(const num_node_t& node);

  // Looks in the tree for unexpanded bits in node. Returns -1 in st of index_t
  // if there was nothing left to expand.
  index_t FindExpansion(int node_idx);

  // Determines if |node_idx| is no longer needed, and removes it and anything else that is now
  // useless.
  void GcNode(int node_idx);
};
}
}
}

#endif  // MEMERNA_SUBOPTIMAL1_NUM_H
