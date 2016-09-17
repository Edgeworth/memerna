#include <algorithm>
#include "fold/suboptimal1_delta.h"

namespace memerna {
namespace fold {
namespace internal {

std::vector<computed_t> Suboptimal1Delta::Run() {
  nodes.push_back({{gext[0][EXT], {0, -1, EXT}}});
  q.push_back(0);
  while (!q.empty()) {
    int node_idx = q.back();
    q.pop_back();
    auto& node = nodes[node_idx];
    assert(node.exp.energy <= max_energy);

    index_t to_expand = node.exp.to_expand;
    // Nothing to do, so either we are done or we need to look up the tree for more nodes to expand.
    if (to_expand.st == -1) {
      assert(node.exp.unexpanded.st == -1);
      to_expand = FindExpansion(node_idx);

      if (to_expand.st == -1) {
        finished.push_back(node_idx);
        continue;
      }
    }

    delta_node_t child = node;
    child.parent = node_idx;
    energy_t base_energy = node.exp.energy;
    const auto& exps = GetExpansions(to_expand);
    for (const auto& exp : exps) {
      child.exp = exp;
      child.exp.energy += base_energy;
      // Since this list is sorted we can break if the energy gets too high.
      if (child.exp.energy > max_energy) break;
      q.push_back(int(nodes.size()));
      nodes.push_back(child);
    }
  }

  auto computeds = ConstructComputedsFromNodes();
  std::sort(computeds.begin(), computeds.end(), computed_energy_comparator_t());
  return computeds;
}

}
}
}
