// Copyright 2023 E.
#ifndef BACKENDS_BASEOPT_SUBOPT_SUBOPT_PERSISTENT_H_
#define BACKENDS_BASEOPT_SUBOPT_SUBOPT_PERSISTENT_H_

#include <algorithm>
#include <cassert>
#include <queue>
#include <utility>
#include <vector>

#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "backends/baseopt/energy/model.h"
#include "backends/baseopt/energy/precomp.h"
#include "backends/common/base/dp.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

// Subopt folding based on a persistent data structure algorithm.
class SuboptPersistent {
 public:
  SuboptPersistent(Primary r, Model::Ptr m, DpState dp, SuboptCfg cfg);

  int Run(const SuboptCallback& fn);

 private:
  struct Node {
    // Index of the parent DfsState in the expand tree.
    int parent_idx{-1};
    // Index of the expansion this DfsState used w.r.t. the parent state's `to_expand`.
    int parent_expand_idx{-1};
    // Index of the next DfsState who's expansion contains an unexpanded DpIndex we need to process.
    int unexpanded_idx{-1};
    // Index of the expansion to use for `unexpanded_idx`'s DfsState.
    int unexpanded_expand_idx{-1};
    // Index of the child expansion of `to_expand` we should process. This gets updated in place
    // (saves time and memory), which is why we need to keep `parent_expand_idx` around as well.
    int expand_idx{0};
    // DpIndex whose child expansions we are processing.
    DpIndex to_expand{};
  };

  Primary r_;
  Model::Ptr m_;
  Precomp pc_;
  DpState dp_;
  SuboptCfg cfg_;
  SuboptResult res_;

  std::vector<std::vector<Expansion>> cache_;
  std::vector<Node> q_;
  std::priority_queue<std::tuple<Energy, int>> pq_;

  std::pair<Energy, int> RunInternal();

  // Computes the suboptimal folding for the given subpath and puts it into `res_`.
  void GenerateResult(int idx);

  const std::vector<Expansion>& GetExpansion(const DpIndex& to_expand) {
    // We request the expansions of an index multiple times when we find the
    // next sibling of a node after coming back up during the DFS.
    auto idx = to_expand.LinearIndex(r_.size());
    if (cache_[idx].empty()) {
      // Need to generate the full way to delta so we can properly set `next_seen`.
      auto exps = GenerateExpansions(to_expand, cfg_.delta);
      std::sort(exps.begin(), exps.end());
      assert(!exps.empty());
      cache_[idx] = std::move(exps);
    }
    return cache_[idx];
  }

  // Generates expansions for the given index, given that the extra energy over the best choice
  // can't be more than `delta`.
  [[nodiscard]] std::vector<Expansion> GenerateExpansions(
      const DpIndex& to_expand, Energy delta) const;
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_SUBOPT_SUBOPT_PERSISTENT_H_
