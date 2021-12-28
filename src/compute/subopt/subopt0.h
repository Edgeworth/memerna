// Copyright 2016 E.
#ifndef COMPUTE_SUBOPT_SUBOPT0_H_
#define COMPUTE_SUBOPT_SUBOPT0_H_

#include <set>
#include <vector>

#include "compute/constants.h"
#include "compute/dp.h"
#include "compute/energy/structure.h"
#include "compute/mfe/globals.h"
#include "compute/subopt/subopt.h"
#include "model/parsing.h"

namespace mrna::subopt::internal {

class Suboptimal0 {
 public:
  Suboptimal0(Energy delta_, int num)
      : max_energy(delta_ == -1 ? CAP_E : mfe::internal::gext[0][EXT] + delta_),
        max_structures(num == -1 ? MAX_STRUCTURES : num) {
    verify(max_structures > 0, "must request at least one structure");
  }
  int Run(SuboptimalCallback fn);

 private:
  struct Node {
    // State should be fully defined by |not_yet_expanded|, |history|, and |base_ctds| which denote
    // what it has done so far, and what it can do from now.
    std::vector<Index> not_yet_expanded;
    std::vector<Index> history;
    std::vector<int16_t> p;
    std::vector<Ctd> base_ctds;
    Energy energy;  // Stores the minimum energy this state could have.

    bool operator<(const Node& o) const { return energy < o.energy; }
  };

  const Energy max_energy;
  const int max_structures;
  // This node is where we build intermediate results to be pushed onto the queue.
  Node curnode;
  std::multiset<Node> finished;
  std::multiset<Node> q;

  void PruneInsert(std::multiset<Node>& prune, const Node& node) {
    if (node.energy <= max_energy) {
      if (static_cast<int>(prune.size()) >= max_structures && (--prune.end())->energy > node.energy)
        prune.erase(--prune.end());
      if (static_cast<int>(prune.size()) < max_structures) prune.insert(node);
    }
  }

  // Creates and inserts a new node with energy |energy| that doesn't
  // need to expand any more ranges than it currently has.
  void Expand(Energy energy) {
    curnode.energy = energy;
    PruneInsert(q, curnode);
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the given ranges.
  void Expand(Energy energy, Index nye) {
    curnode.not_yet_expanded.push_back(nye);
    curnode.energy = energy;
    PruneInsert(q, curnode);
    curnode.not_yet_expanded.pop_back();
  }

  void Expand(Energy energy, Index nye, IndexCtd ctd_idx) {
    curnode.base_ctds[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye);
    curnode.base_ctds[ctd_idx.idx] = CTD_NA;
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the two given ranges.
  void Expand(Energy energy, Index nye0, Index nye1) {
    curnode.not_yet_expanded.push_back(nye0);
    curnode.not_yet_expanded.push_back(nye1);
    curnode.energy = energy;
    PruneInsert(q, curnode);
    curnode.not_yet_expanded.pop_back();
    curnode.not_yet_expanded.pop_back();
  }

  void Expand(Energy energy, Index nye0, Index nye1, IndexCtd ctd_idx) {
    curnode.base_ctds[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye0, nye1);
    curnode.base_ctds[ctd_idx.idx] = CTD_NA;
  }

  void Expand(Energy energy, Index nye0, Index nye1, IndexCtd ctd_idx0, IndexCtd ctd_idx1) {
    curnode.base_ctds[ctd_idx0.idx] = ctd_idx0.ctd;
    curnode.base_ctds[ctd_idx1.idx] = ctd_idx1.ctd;
    Expand(energy, nye0, nye1);
    curnode.base_ctds[ctd_idx0.idx] = CTD_NA;
    curnode.base_ctds[ctd_idx1.idx] = CTD_NA;
  }
};

}  // namespace mrna::subopt::internal

#endif  // COMPUTE_SUBOPT_SUBOPT0_H_
