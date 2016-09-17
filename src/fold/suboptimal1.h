#ifndef MEMERNA_SUBOPTIMAL1_H
#define MEMERNA_SUBOPTIMAL1_H

#include <set>
#include <algorithm>
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
        q([this](int a, int b) {return NodeComparator(a, b);}) {
    verify_expr(max_structures > 0, "must request at least one structure");
  }
  std::vector<computed_t> Run();

private:
  struct expand_t {
    expand_t() = delete;
    expand_t(energy_t energy_) : energy(energy_) {}
    expand_t(energy_t energy_, const index_t& to_expand_)
        : energy(energy_), to_expand(to_expand_) {}
    expand_t(energy_t energy_, const index_t& to_expand_, const index_t& unexpanded_)
        : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_) {}
    expand_t(energy_t energy_, const index_t& to_expand_, const ctd_idx_t& ctd0_)
        : energy(energy_), to_expand(to_expand_), ctd0(ctd0_) {}
    expand_t(energy_t energy_, const index_t& to_expand_,
        const index_t& unexpanded_, const ctd_idx_t& ctd0_)
        : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_), ctd0(ctd0_) {}
    expand_t(energy_t energy_, const index_t& to_expand_,
        const index_t& unexpanded_, const ctd_idx_t& ctd0_, const ctd_idx_t& ctd1_)
        : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_), ctd0(ctd0_), ctd1(ctd1_) {}
    energy_t energy;
    index_t to_expand, unexpanded;  // st is -1 if this does not exist
    ctd_idx_t ctd0, ctd1;

    bool operator<(const expand_t& o) const {return energy < o.energy;}
  };

  struct node_t {
    node_t() = delete;
    node_t(const expand_t& exp_)
        : exp(exp_), child_count(0), parent(-1),
          exp_st(0), exp_en(-1), cur_anc(0) {}
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
    int child_count, parent, exp_st, exp_en, cur_anc;
  };

  const energy_t max_energy;
  const int max_structures;
  std::vector<node_t> nodes;
  // This node is where we build intermediate results to be pushed onto the queue.
  std::vector<int> finished;
  std::vector<int> free_space;
  std::unordered_map<index_t, std::vector<expand_t>> cache;
  std::set<int, std::function<bool(int, int)>> q;

  bool NodeComparator(int a, int b) const {
    if (nodes[a].exp.energy != nodes[b].exp.energy)
      return nodes[a].exp.energy < nodes[b].exp.energy;
    if (a != b) return a < b;
    return false;
  }

  bool InsertQ(const node_t& node) {
    bool inserted = false;
    if (node.exp.energy <= max_energy) {
      int gc_node = -1;
      if (int(q.size()) >= max_structures && nodes[*(--q.end())].exp.energy > node.exp.energy) {
        gc_node = *(--q.end());
        q.erase(--q.end());
      }
      if (int(q.size()) < max_structures) {
        if (!free_space.empty()) {
          int node_idx = free_space.back();
          free_space.pop_back();
          nodes[node_idx] = node;
          q.insert(node_idx);
        } else {
          nodes.push_back(node);
          q.insert(int(nodes.size() - 1));
        }
        // Update child counts.
        if (node.parent != -1)
          ++(nodes[node.parent].child_count);
        inserted = true;
      }
      // We have to GC after doing the insert because of this one weird case:
      // We inserted a child previously whose parent has just that as its one child.
      // We then removed that child so we can insert |node|, whose parent is the same.
      // If we had GC'd before, we would have removed the common parent.
      if (gc_node != -1) GcNode(gc_node);
    }
    return inserted;
  }

  // Looks in the tree for unexpanded bits in node. Returns -1 in st of index_t
  // if there was nothing left to expand.
  index_t FindExpansion(int node_idx);

  // Traverses up the tree and reconstructs a computed_t.
  computed_t ReconstructComputed(int node_idx);

  // Determines if |node_idx| is no longer needed, and removes it and anything else that is now useless.
  void GcNode(int node_idx);

  std::vector<expand_t> GenerateExpansions(const index_t& to_expand);
};

}
}
}

#endif   // MEMERNA_SUBOPTIMAL1_H
