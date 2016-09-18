#ifndef MEMERNA_SUBOPTIMAL1_BASE_H
#define MEMERNA_SUBOPTIMAL1_BASE_H

#include <algorithm>
#include "common.h"
#include "fold/fold_internal.h"

namespace memerna {
namespace fold {
namespace internal {

// TODO try to make this smaller?
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

std::vector<expand_t> GenerateExpansions(const index_t& to_expand);

template<typename Node>
class Suboptimal1Base {
protected:
  std::vector<Node> nodes;
  std::vector<int> finished;  // TODO get rid of this, don't need it with new reconstruction method.

  // Looks in the tree for unexpanded bits in node. Returns -1 in st of index_t
  // if there was nothing left to expand.
  index_t FindExpansion(int node_idx) {
    auto& node = nodes[node_idx];
    while (1) {
      if (node.exp_en == node_idx) {
        // Finished state - when |exp_en| is us, then we have been expanded, and everything above us.
        return {};
      } else {
        // Need to find the next node to expand in the tree. Keep looking up the tree while
        // we haven't found a node which we can expand, and while we are still in range.
        while (node.cur_anc != node.exp_en && nodes[node.cur_anc].exp.unexpanded.st == -1) {
          node.cur_anc = nodes[node.cur_anc].parent;
          assert(node.cur_anc == -1 || nodes[node.cur_anc].child_count != 0);
        }
        if (node.cur_anc == node.exp_en) {
          // Case: We did not find a node. In this case, update expand range to be from us to exp_st.
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

  std::vector<computed_t> ConstructComputedsFromNodes() {
    std::vector<computed_t> ret;
    for (int node_idx : finished)
      ret.push_back(ReconstructComputed(node_idx));
    return ret;
  }

  const std::vector<expand_t>& GetExpansions(const index_t& to_expand) {
    auto iter = cache.find(to_expand);
    if (iter == cache.end()) {
      auto exps = GenerateExpansions(to_expand);
      std::sort(exps.begin(), exps.end());
      auto res = cache.emplace(to_expand, std::move(exps));
      assert(res.second && res.first != cache.end());
      iter = res.first;
    }
    return iter->second;
  }

private:
  // This node is where we build intermediate results to be pushed onto the queue.
  std::unordered_map<index_t, std::vector<expand_t>> cache;

  computed_t ReconstructComputed(int node_idx) {
    assert(nodes[node_idx].exp.to_expand.st == -1 && nodes[node_idx].exp_en == node_idx);
    computed_t computed(gr);
    computed.energy = nodes[node_idx].exp.energy;
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
      if (node.exp.ctd0.idx != -1)
        computed.base_ctds[node.exp.ctd0.idx] = node.exp.ctd0.ctd;
      if (node.exp.ctd1.idx != -1)
        computed.base_ctds[node.exp.ctd1.idx] = node.exp.ctd1.ctd;
      node_idx = node.parent;
    }
    return computed;
  }
};

}
}
}

#endif  // MEMERNA_SUBOPTIMAL1_BASE_H
