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
  Suboptimal0(Primary r, energy::EnergyModel em, Energy delta, int num);

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

  Primary r_;
  energy::EnergyModel em_;

  const Energy max_energy_;
  const int max_structures;
  // This node is where we build intermediate results to be pushed onto the queue.
  Node curnode_;
  std::multiset<Node> finished_;
  std::multiset<Node> q_;

  void PruneInsert(std::multiset<Node>& prune, const Node& node) {
    if (node.energy <= max_energy_) {
      if (static_cast<int>(prune.size()) >= max_structures && (--prune.end())->energy > node.energy)
        prune.erase(--prune.end());
      if (static_cast<int>(prune.size()) < max_structures) prune.insert(node);
    }
  }

  // Creates and inserts a new node with energy |energy| that doesn't
  // need to expand any more ranges than it currently has.
  void Expand(Energy energy) {
    curnode_.energy = energy;
    PruneInsert(q_, curnode_);
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the given ranges.
  void Expand(Energy energy, Index nye) {
    curnode_.not_yet_expanded.push_back(nye);
    curnode_.energy = energy;
    PruneInsert(q_, curnode_);
    curnode_.not_yet_expanded.pop_back();
  }

  void Expand(Energy energy, Index nye, IndexCtd ctd_idx) {
    curnode_.base_ctds[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye);
    curnode_.base_ctds[ctd_idx.idx] = CTD_NA;
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the two given ranges.
  void Expand(Energy energy, Index nye0, Index nye1) {
    curnode_.not_yet_expanded.push_back(nye0);
    curnode_.not_yet_expanded.push_back(nye1);
    curnode_.energy = energy;
    PruneInsert(q_, curnode_);
    curnode_.not_yet_expanded.pop_back();
    curnode_.not_yet_expanded.pop_back();
  }

  void Expand(Energy energy, Index nye0, Index nye1, IndexCtd ctd_idx) {
    curnode_.base_ctds[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye0, nye1);
    curnode_.base_ctds[ctd_idx.idx] = CTD_NA;
  }

  void Expand(Energy energy, Index nye0, Index nye1, IndexCtd ctd_idx0, IndexCtd ctd_idx1) {
    curnode_.base_ctds[ctd_idx0.idx] = ctd_idx0.ctd;
    curnode_.base_ctds[ctd_idx1.idx] = ctd_idx1.ctd;
    Expand(energy, nye0, nye1);
    curnode_.base_ctds[ctd_idx0.idx] = CTD_NA;
    curnode_.base_ctds[ctd_idx1.idx] = CTD_NA;
  }
};

}  // namespace mrna::subopt::internal

#endif  // COMPUTE_SUBOPT_SUBOPT0_H_
