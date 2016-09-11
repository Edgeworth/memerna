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
  Suboptimal0(energy_t max_energy_, int max_structures_);
  std::vector<computed_t> Run();

private:
  struct suboptimal_node_t {
    // TODO Since nodes form a tree, can save memory on these two vectors.
    // TODO Will also save time in the comparison.
    // State should be fully defined by |not_yet_expanded| and |history|, which denote
    // what it has done so far, and what it can do from now.
    std::vector<index_t> not_yet_expanded;
    std::vector<index_t> history;
    std::vector<int> p;
    std::vector<Ctd> base_ctds;
    energy_t energy;  // Stores the minimum energy this state could have.

    bool operator<(const suboptimal_node_t& o) const {
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
  suboptimal_node_t curnode;
  std::set<suboptimal_node_t> finished;
  std::set<suboptimal_node_t> q;

  void PrintNodeDebug(const suboptimal_node_t& node) {
    printf("Node - not yet expanded:\n  ");
    for (const auto& idx : node.not_yet_expanded)
      printf("(%d, %d, %d), ", idx.st, idx.en, idx.a);
    printf("\nhistory:\n  ");
    for (const auto& idx : node.history)
      printf("(%d, %d, %d), ", idx.st, idx.en, idx.a);
    printf("\npairs so far: %s\nbase ctds:\n  ", parsing::PairsToDotBracket(node.p).c_str());
    for (const auto& ctd : node.base_ctds) {
      printf("%s, ", energy::CtdToName(ctd));
    }
    printf("\nenergy: %d\n", node.energy);
  }

  void PruneInsert(std::set<suboptimal_node_t>& prune, const suboptimal_node_t& node) {
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

  void Expand(energy_t energy, index_t nye, std::pair<Ctd, int> ctd_idx) {
    curnode.base_ctds[ctd_idx.second] = ctd_idx.first;
    Expand(energy, nye);
    curnode.base_ctds[ctd_idx.second] = CTD_NA;
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

  void Expand(energy_t energy, index_t nye0, index_t nye1, std::pair<Ctd, int> ctd_idx) {
    curnode.base_ctds[ctd_idx.second] = ctd_idx.first;
    Expand(energy, nye0, nye1);
    curnode.base_ctds[ctd_idx.second] = CTD_NA;
  }

  void Expand(energy_t energy, index_t nye0, index_t nye1,
      std::pair<Ctd, int> ctd_idx0, std::pair<Ctd, int> ctd_idx1) {
    curnode.base_ctds[ctd_idx0.second] = ctd_idx0.first;
    curnode.base_ctds[ctd_idx1.second] = ctd_idx1.first;
    Expand(energy, nye0, nye1);
    curnode.base_ctds[ctd_idx0.second] = CTD_NA;
    curnode.base_ctds[ctd_idx1.second] = CTD_NA;
  }
};

}
}
}

#endif   // MEMERNA_SUBOPTIMAL0_H