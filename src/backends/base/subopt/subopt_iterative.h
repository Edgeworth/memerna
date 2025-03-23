// Copyright 2016 Eliot Courtney.
#ifndef BACKENDS_BASE_SUBOPT_SUBOPT_ITERATIVE_H_
#define BACKENDS_BASE_SUBOPT_SUBOPT_ITERATIVE_H_

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/base/energy/precomp.h"
#include "backends/common/base/dp.h"
#include "backends/common/expansion_cache.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

template <bool UseLru>
class SuboptIterative {
 public:
  SuboptIterative(Primary r, Model::Ptr m, DpState dp, SuboptCfg cfg);

  int Run(const SuboptCallback& fn);

 private:
  struct Node {
    // Index of the child expansion of `to_expand` we should process.
    int expand_idx = {0};
    // DpIndex whose child expansions we are processing. -1 means empty
    DpIndex to_expand{};
    // Stores whether this node's `to_expand` was from `unexpanded_` and needs
    // to be replaced when going back up the DFS stack.
    bool should_unexpand = {false};
  };

  Primary r_;
  Model::Ptr m_;
  Precomp pc_;
  DpState dp_;
  SuboptCfg cfg_;

  ExpansionCache<DpIndex, Expansion, UseLru> cache_;
  std::vector<Node> q_;

  // Incremental state. Holds the current partial structure.
  SuboptResult res_;
  // Incremental state - holds unexpanded Indexes for the current partial structure.
  std::vector<DpIndex> unexpanded_;

  std::pair<int, Energy> RunInternal(
      const SuboptCallback& fn, Energy delta, bool exact_energy, int max);

  const std::vector<Expansion>& GetExpansion(const DpIndex& to_expand) {
    // We request the expansions of an index multiple times when we find the
    // next sibling of a node after coming back up during the DFS.
    auto key = to_expand.LinearIndex(r_.size());
    if (const auto& val = cache_.Get(key); !val.empty()) return val;

    // Need to generate the full way to delta so we can properly set `next_seen`.
    auto exps = GenerateExpansions(to_expand, cfg_.delta);
    std::sort(exps.begin(), exps.end());
    assert(!exps.empty());
    return cache_.Insert(key, std::move(exps));
  }

  // Generates expansions for the given index, given that the extra energy over the best choice
  // can't be more than `delta`.
  [[nodiscard]] std::vector<Expansion> GenerateExpansions(
      const DpIndex& to_expand, Energy delta) const;
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_SUBOPT_SUBOPT_ITERATIVE_H_
