#include "fold/suboptimal1_num.h"

namespace memerna {
namespace fold {
namespace internal {

bool Suboptimal1Num::InsertQ(const num_node_t& node) {
  bool inserted = false;
  if (node.exp.energy <= delta) {
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

void Suboptimal1Num::GcNode(int node_idx) {
  while (nodes[node_idx].child_count == 0) {
    free_space.push_back(node_idx);
    node_idx = nodes[node_idx].parent;
    if (node_idx == -1) break;
    --(nodes[node_idx].child_count);
  }
}

std::vector<computed_t> Suboptimal1Num::Run() {
  InsertQ({{0, {0, -1, EXT}}});
  while (!q.empty()) {
    int node_idx = *q.begin();
    auto& node = nodes[node_idx];
    q.erase(q.begin());
    assert(node.exp.energy <= delta);

    index_t to_expand = node.exp.to_expand;
    // Nothing to do, so either we are done or we need to look up the tree for more nodes to expand.
    if (to_expand.st == -1) {
      assert(node.exp.unexpanded.st == -1);
      to_expand = FindExpansion(node_idx);

      if (to_expand.st == -1) {
        finished.push_back(node_idx);
        // Since we add into |finished| in order, if we have enough structures, exit.
        if (max_structures == int(finished.size())) break;
        else continue;
      }
    }

    // Make a copy now, since the parent might get GC'd.
    num_node_t child = node;
    child.parent = node_idx;
    energy_t base_energy = node.exp.energy;
    const auto& exps = GetExpansions(to_expand);
    for (const auto& exp : exps) {
      child.exp = exp;
      child.exp.energy += base_energy;
      // Since this list is sorted we can break if the energy gets too high.
      if (!InsertQ(child))
        break;
    }
  }

  // Should be already sorted.
  return ConstructComputedsFromNodes();
}

}
}
}
