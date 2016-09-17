#ifndef MEMERNA_SUBOPTIMAL1_NUM_H
#define MEMERNA_SUBOPTIMAL1_NUM_H

#include <set>
#include "common.h"
#include "fold/suboptimal1_base.h"
#include "fold/fold_internal.h"

namespace memerna {
namespace fold {
namespace internal {

struct num_node_t {
  num_node_t() = delete;
  num_node_t(const expand_t& exp_)
      : exp(exp_), child_count(0), parent(-1),
        exp_st(0), exp_en(-1), cur_anc(0) {}
  expand_t exp;
  // See Suboptimal1Delta for explanation.
  int child_count, parent, exp_st, exp_en, cur_anc;
};

class Suboptimal1Num : public Suboptimal1Base<num_node_t> {
public:
  Suboptimal1Num(energy_t delta_, int num)
      : max_energy(delta_ == -1 ? CAP_E : delta_ + gext[0][EXT]),
        max_structures(num == -1 ? std::numeric_limits<int>::max() / 4 : num),
        q([this](int a, int b) {return NodeComparator(a, b);}) {
    verify_expr(max_structures > 0, "must request at least one structure");
  }

  std::vector<computed_t> Run();

private:
  const energy_t max_energy;
  const int max_structures;
  std::vector<int> free_space;  // Free space from stuff that has been GC'd.
  std::set<int, std::function<bool(int, int)>> q;  // Queue of nodes ready to be expanded.

  bool NodeComparator(int a, int b) const {
    if (nodes[a].exp.energy != nodes[b].exp.energy)
      return nodes[a].exp.energy < nodes[b].exp.energy;
    return a < b;
  }

  bool InsertQ(const num_node_t& node);

  // Determines if |node_idx| is no longer needed, and removes it and anything else that is now useless.
  void GcNode(int node_idx);
};


}
}
}

#endif  // MEMERNA_SUBOPTIMAL1_NUM_H
