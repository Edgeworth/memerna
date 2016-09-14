#ifndef MEMERNA_SUBOPTIMAL1_H
#define MEMERNA_SUBOPTIMAL1_H

#include <set>
#include "common.h"
#include "parsing.h"
#include "fold/fold_internal.h"
#include "energy/structure.h"

namespace memerna {
namespace fold {
namespace internal {

class Suboptimal1 {
public:
  Suboptimal1(energy_t max_energy_, int max_structures_)
    : max_energy(max_energy_), max_structures(max_structures_),
      finished([this](int a, int b) {return NodeComparator(a, b);}),
      q([this](int a, int b) {return NodeComparator(a, b);}) {
    verify_expr(max_structures > 0, "must request at least one structure");
  }
  std::vector<computed_t> Run();

private:
  struct node_t {
    // TODO try to reduce size of this? - bitfields, etc - 28 or so bytes might be possible
    // Size limited heap
    // TODO Cache results of expanding - use GC to keep track of unexpanded's that will
    // be expanded multiple times?

    energy_t energy;
    index_t to_expand;  // st is -1 if this does not exist
    index_t unexpanded;  // st is -1 if this does not exist
    int child_count;  // For GC.
    int parent;
    // Method for expanding ancestors |unexpanded|s:
    // Keep track of lowest ancestor which it and everything above is expanded, |expand_en|.
    // Cur ancestor starts at |expand_st| tracks the next ancestor to expand, and goes up the tree.
    // Once it reaches |expand_en|, update |expand_en| to |expand_st|, and
    // |expand_st| and |cur_ancestor|to the current node.
    // [expand_st, expand_en) is the half-open range saying which nodes we should look
    // for things to expand in.
    int expand_st;
    int expand_en;
    int cur_ancestor;
    std::pair<Ctd, int> ctd0, ctd1;
  };

  const energy_t max_energy;
  const int max_structures;
  // This node is where we build intermediate results to be pushed onto the queue.
  std::set<int, std::function<bool(int, int)>> finished;
  std::set<int, std::function<bool(int, int)>> q;
  std::vector<int> free_space;
  std::vector<node_t> nodes;

  bool NodeComparator(int a, int b) const {
    if (nodes[a].energy != nodes[b].energy) return nodes[a].energy < nodes[b].energy;
    if (a != b) return a < b;
    return false;
  }

  void PrintNodeDebug(int node_idx) {
    const auto& node = nodes[node_idx];
    printf("Node %d, energy %d, parent: %d, expand_st: %d, expand_en: %d, cur_anc: %d:\n",
        node_idx, node.energy, node.parent, node.expand_st, node.expand_en, node.cur_ancestor);
    printf("  %d %d %d; %d %d %d\n", node.to_expand.st, node.to_expand.en, node.to_expand.a,
        node.unexpanded.st, node.unexpanded.en, node.unexpanded.a);
  }

  void InsertFinished(int node_idx) {
    const auto& node = nodes[node_idx];
    if (node.energy <= max_energy) {
      if (int(finished.size()) >= max_structures && nodes[*(--finished.end())].energy > node.energy) {
        GcNode(*(--finished.end()));
        finished.erase(--finished.end());
      }
      if (int(finished.size()) < max_structures)
        finished.insert(node_idx);
    }
  }

  void InsertQ(const node_t& node) {
    if (node.energy <= max_energy) {
      int gc_node = -1;
      if (int(q.size()) >= max_structures && nodes[*(--q.end())].energy > node.energy) {
        gc_node = *(--q.end());
        q.erase(--q.end());
      }
      if (int(q.size()) < max_structures) {
        if (!free_space.empty()) {
          int node_idx = free_space.back();
          free_space.pop_back();
          assert(q.count(node_idx) == 0);
          nodes[node_idx] = node;
          q.insert(node_idx);
        } else {
          nodes.push_back(node);
          q.insert(int(nodes.size() - 1));
        }
        // Update child counts.
        if (node.parent != -1)
          ++(nodes[node.parent].child_count);
      }
      // We have to GC after doing the insert because of this one weird case:
      // We inserted a child previously whose parent has just that as its one child.
      // We then removed that child so we can insert |node|, whose parent is the same.
      // If we had GC'd before, we would have removed the common parent.
      if (gc_node != -1) GcNode(gc_node);
    }
  }

  // Creates and inserts a new node with energy |energy| that doesn't
  // need to expand any more ranges than it currently has.
  void Expand(int parent, energy_t energy) {
    node_t node = {energy, {-1, -1, -1}, {-1, -1, -1}, 0, parent,
        nodes[parent].expand_st, nodes[parent].expand_en, nodes[parent].cur_ancestor,
        {CTD_NA, -1}, {CTD_NA, -1}
    };
    InsertQ(node);
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the given ranges.
  void Expand(int parent, energy_t energy, index_t nye) {
    node_t node = {energy, nye, {-1, -1, -1}, 0, parent,
        nodes[parent].expand_st, nodes[parent].expand_en, nodes[parent].cur_ancestor,
        {CTD_NA, -1}, {CTD_NA, -1}
    };
    InsertQ(node);
  }

  void Expand(int parent, energy_t energy, index_t nye, std::pair<Ctd, int> ctd_idx) {
    node_t node = {energy, nye, {-1, -1, -1}, 0, parent,
        nodes[parent].expand_st, nodes[parent].expand_en, nodes[parent].cur_ancestor,
        ctd_idx, {CTD_NA, -1}
    };
    InsertQ(node);
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the two given ranges.
  void Expand(int parent, energy_t energy, index_t nye0, index_t nye1) {
    node_t node = {energy, nye0, nye1, 0, parent,
        nodes[parent].expand_st, nodes[parent].expand_en, nodes[parent].cur_ancestor,
        {CTD_NA, -1}, {CTD_NA, -1}
    };
    InsertQ(node);
  }

  void Expand(int parent, energy_t energy, index_t nye0, index_t nye1, std::pair<Ctd, int> ctd_idx) {
    node_t node = {energy, nye0, nye1, 0, parent,
        nodes[parent].expand_st, nodes[parent].expand_en, nodes[parent].cur_ancestor,
        ctd_idx, {CTD_NA, -1}
    };
    InsertQ(node);
  }

  void Expand(int parent, energy_t energy, index_t nye0, index_t nye1,
      std::pair<Ctd, int> ctd_idx0, std::pair<Ctd, int> ctd_idx1) {
    node_t node = {energy, nye0, nye1, 0, parent,
        nodes[parent].expand_st, nodes[parent].expand_en, nodes[parent].cur_ancestor,
        ctd_idx0, ctd_idx1
    };
    InsertQ(node);
  }

  // Looks in the tree for unexpanded bits in node. Returns -1 in st of index_t
  // if there was nothing left to expand.
  index_t FindExpansion(int node_idx);

  // Traverses up the tree and reconstructs a computed_t.
  computed_t ReconstructComputed(int node_idx);

  // Determines if |node_idx| is no longer needed, and removes it and anything else that is now useless.
  void GcNode(int node_idx);
};

}
}
}

#endif   // MEMERNA_SUBOPTIMAL1_H
