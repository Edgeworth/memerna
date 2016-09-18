#ifndef MEMERNA_SUBOPTIMAL1_DELTA_H
#define MEMERNA_SUBOPTIMAL1_DELTA_H

#include "common.h"
#include "fold/suboptimal1_base.h"
#include "fold/fold_internal.h"

namespace memerna {
namespace fold {
namespace internal {

struct delta_node_t {
  delta_node_t() = delete;
  delta_node_t(const expand_t& exp_)
      : exp(exp_), parent(-1), exp_st(0), exp_en(-1), cur_anc(0) {}
  // TODO try to reduce size of this? - bitfields, etc - 28 or so bytes might be possible
  expand_t exp;
  // Book-keeping vars:
  // Method for expanding ancestors |unexpanded|s:
  // Keep track of lowest ancestor which it and everything above is expanded, |exp_en|.
  // Cur ancestor starts at |exp_st| tracks the next ancestor to expand, and goes up the tree.
  // Once it reaches |exp_en|, update |exp_en| to |exp_st|, and
  // |exp_st| and |cur_anc|to the current node.
  // [exp_st, exp_en) is the half-open range saying which nodes we should look
  // for things to expand in.
  int parent, exp_st, exp_en, cur_anc;
};

class Suboptimal1Delta : public Suboptimal1Base<delta_node_t> {
public:
  Suboptimal1Delta(energy_t delta_) : Suboptimal1Base(delta_) {
    verify_expr(delta_ >= 0, "energy delta must be non-negative");
  }

  std::vector<computed_t> Run();

private:
  std::vector<int> q;  // Queue of nodes ready to be expanded.
};

}
}
}

#endif  // MEMERNA_SUBOPTIMAL1_DELTA_H
