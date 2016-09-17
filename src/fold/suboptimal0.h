#ifndef MEMERNA_SUBOPTIMAL0_H
#define MEMERNA_SUBOPTIMAL0_H

#include <set>
#include "common.h"
#include "parsing.h"
#include "fold/fold_internal.h"
#include "energy/structure.h"

namespace memerna {
namespace fold {
namespace internal {

class Suboptimal0 {
public:
  Suboptimal0(energy_t max_energy_, int max_structures_)
      : max_energy(max_energy_), max_structures(max_structures_) {
    verify_expr(max_structures > 0, "must request at least one structure");
  }
  std::vector<computed_t> Run();

private:
  struct node_t {
    // State should be fully defined by |not_yet_expanded|, |history|, and |base_ctds| which denote
    // what it has done so far, and what it can do from now.
    std::vector<index_t> not_yet_expanded;
    std::vector<index_t> history;
    std::vector<int> p;
    std::vector<Ctd> base_ctds;
    energy_t energy;  // Stores the minimum energy this state could have.

    bool operator<(const node_t& o) const {
      if (energy != o.energy) return energy < o.energy;
      if (not_yet_expanded != o.not_yet_expanded) return not_yet_expanded < o.not_yet_expanded;
      if (history != o.history) return history < o.history;
      if (base_ctds != o.base_ctds) return base_ctds < o.base_ctds;
      assert(false);  // Should never happen.
      return false;
    }
  };

  const energy_t max_energy;
  const int max_structures;
  // This node is where we build intermediate results to be pushed onto the queue.
  node_t curnode;
  std::set<node_t> finished;
  std::set<node_t> q;

  void PruneInsert(std::set<node_t>& prune, const node_t& node) {
    if (node.energy <= max_energy) {
      if (int(prune.size()) >= max_structures && (--prune.end())->energy > node.energy)
        prune.erase(--prune.end());
      if (int(prune.size()) < max_structures)
        prune.insert(node);
    }
  }

  // Creates and inserts a new node with energy |energy| that doesn't
  // need to expand any more ranges than it currently has.
  void Expand(energy_t energy) {
    curnode.energy = energy;
    PruneInsert(q, curnode);
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the given ranges.
  void Expand(energy_t energy, index_t nye) {
    curnode.not_yet_expanded.push_back(nye);
    curnode.energy = energy;
    PruneInsert(q, curnode);
    curnode.not_yet_expanded.pop_back();
  }

  void Expand(energy_t energy, index_t nye, ctd_idx_t ctd_idx) {
    curnode.base_ctds[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye);
    curnode.base_ctds[ctd_idx.idx] = CTD_NA;
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the two given ranges.
  void Expand(energy_t energy, index_t nye0, index_t nye1) {
    curnode.not_yet_expanded.push_back(nye0);
    curnode.not_yet_expanded.push_back(nye1);
    curnode.energy = energy;
    PruneInsert(q, curnode);
    curnode.not_yet_expanded.pop_back();
    curnode.not_yet_expanded.pop_back();
  }

  void Expand(energy_t energy, index_t nye0, index_t nye1, ctd_idx_t ctd_idx) {
    curnode.base_ctds[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye0, nye1);
    curnode.base_ctds[ctd_idx.idx] = CTD_NA;
  }

  void Expand(energy_t energy, index_t nye0, index_t nye1,
      ctd_idx_t ctd_idx0, ctd_idx_t ctd_idx1) {
    curnode.base_ctds[ctd_idx0.idx] = ctd_idx0.ctd;
    curnode.base_ctds[ctd_idx1.idx] = ctd_idx1.ctd;
    Expand(energy, nye0, nye1);
    curnode.base_ctds[ctd_idx0.idx] = CTD_NA;
    curnode.base_ctds[ctd_idx1.idx] = CTD_NA;
  }
};

}
}
}

#endif   // MEMERNA_SUBOPTIMAL0_H
