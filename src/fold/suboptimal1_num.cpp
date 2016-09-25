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
      if (node.parent != -1) ++(nodes[node.parent].child_count);
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

void Suboptimal1Num::Run(std::function<void(const computed_t&)> fn) {
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
        if (max_structures == int(finished.size()))
          break;
        else
          continue;
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
      if (!InsertQ(child)) break;
    }
  }

  // Should be already sorted
  // TODO depends on stuff but at least put this fn call inside construct
  // TODO computeds and don't use so much memory
  for (const auto& computed : ConstructComputedsFromNodes())
    fn(computed);
}

computed_t Suboptimal1Num::ReconstructComputed(int node_idx) {
  assert(nodes[node_idx].exp.to_expand.st == -1 && nodes[node_idx].exp_en == node_idx);
  computed_t computed(gr);
  computed.energy = nodes[node_idx].exp.energy + gext[0][EXT];
  while (node_idx != -1) {
    const auto& node = nodes[node_idx];
    if (node.exp.to_expand.en != -1 && node.exp.to_expand.a == DP_P) {
      computed.s.p[node.exp.to_expand.st] = node.exp.to_expand.en;
      computed.s.p[node.exp.to_expand.en] = node.exp.to_expand.st;
    }
    // Also need to grab unexpanded ones - these will have been expanded, but
    // they themselves don't get replicated as to_expands in the tree.
    if (node.exp.unexpanded.en != -1 && node.exp.unexpanded.a == DP_P) {
      computed.s.p[node.exp.unexpanded.st] = node.exp.unexpanded.en;
      computed.s.p[node.exp.unexpanded.en] = node.exp.unexpanded.st;
    }
    if (node.exp.ctd0.idx != -1) computed.base_ctds[node.exp.ctd0.idx] = node.exp.ctd0.ctd;
    if (node.exp.ctd1.idx != -1) computed.base_ctds[node.exp.ctd1.idx] = node.exp.ctd1.ctd;
    node_idx = node.parent;
  }
  return computed;
}

index_t Suboptimal1Num::FindExpansion(int node_idx) {
  auto& node = nodes[node_idx];
  while (1) {
    if (node.exp_en == node_idx) {
      // Finished state - when |exp_en| is us, then we have been expanded, and everything above
      // us.
      return {};
    } else {
      // Need to find the next node to expand in the tree. Keep looking up the tree while
      // we haven't found a node which we can expand, and while we are still in range.
      while (node.cur_anc != node.exp_en && nodes[node.cur_anc].exp.unexpanded.st == -1)
        node.cur_anc = nodes[node.cur_anc].parent;
      if (node.cur_anc == node.exp_en) {
        // Case: We did not find a node. In this case, update expand range to be from us to
        // exp_st.
        node.exp_en = node.exp_st;
        node.exp_st = node_idx;
        node.cur_anc = node_idx;
        // Need to repeat now, and explore the new range.
      } else {
        // Case: We found a node to expand
        auto to_expand = nodes[node.cur_anc].exp.unexpanded;
        node.cur_anc = nodes[node.cur_anc].parent;  // Just used this one, so go up the tree.
        return to_expand;
      }
    }
  }
}

std::vector<computed_t> Suboptimal1Num::ConstructComputedsFromNodes() {
  std::vector<computed_t> ret;
  for (int node_idx : finished)
    ret.push_back(ReconstructComputed(node_idx));
  return ret;
}

}
}
}
