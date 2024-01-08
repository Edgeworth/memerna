// Copyright 2023 Eliot Courtney.
#ifndef MODELS_T04_SUBOPT_SUBOPT_PERSISTENT_H_
#define MODELS_T04_SUBOPT_SUBOPT_PERSISTENT_H_

#include <algorithm>
#include <cassert>
#include <compare>
#include <functional>
#include <queue>
#include <utility>
#include <vector>

#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "models/t04/energy/model.h"
#include "models/t04/energy/precomp.h"
#include "models/t04/mfe/dp.h"
#include "models/t04/trace/trace.h"
#include "util/splaymap.h"

namespace mrna::md::t04 {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

// Suboptimal folding based on a persistent data structure algorithm.
class SuboptPersistent {
 public:
  SuboptPersistent(Primary r, Model::Ptr em, DpState dp, SuboptCfg cfg);

  int Run(const SuboptCallback& fn);

 private:
  struct Node {
    // Index of the parent DfsState in the expand tree.
    int parent_idx = {-1};
    // Index of the expansion this DfsState used w.r.t. the parent state's `to_expand`.
    int parent_expand_idx = {-1};
    // Index of the next DfsState who's expansion contains an unexpanded DpIndex we need to process.
    int unexpanded_idx = {-1};
    // Index of the expansion to use for `unexpanded_idx`'s DfsState.
    int unexpanded_expand_idx = {-1};
    // Index of the child expansion of `to_expand` we should process. This gets updated in place
    // (saves time and memory), which is why we need to keep `parent_expand_idx` around as well.
    int expand_idx = {0};
    // DpIndex whose child expansions we are processing.
    DpIndex to_expand{};
  };

  Primary r_;
  Model::Ptr em_;
  Precomp pc_;
  DpState dp_;
  SuboptCfg cfg_;
  SuboptResult res_;

  SplayMap<DpIndex, std::vector<Expansion>> cache_;
  std::vector<Node> q_;
  std::priority_queue<std::pair<Energy, int>> pq_;

  std::pair<Energy, int> RunInternal();

  // Computes the suboptimal folding for the given subpath and puts it into `res_`.
  void GenerateResult(int idx);

  const std::vector<Expansion>& GetExpansion(const DpIndex& to_expand) {
    // We request the expansions of an index multiple times when we find the
    // next sibling of a node after coming back up during the DFS.
    if (!cache_.Find(to_expand)) {
      // Need to generate the full way to delta so we can properly set `next_seen`.
      auto exps = GenerateExpansions(to_expand, cfg_.delta);
      std::sort(exps.begin(), exps.end());
      [[maybe_unused]] auto res = cache_.Insert(to_expand, std::move(exps));
      assert(res);
    }
    return cache_.Get();
  }

  // Generates expansions for the given index, given that the extra energy over the best choice
  // can't be more than `delta`.
  [[nodiscard]] std::vector<Expansion> GenerateExpansions(
      const DpIndex& to_expand, Energy delta) const;
};

}  // namespace mrna::md::t04

#endif  // MODELS_T04_SUBOPT_SUBOPT_PERSISTENT_H_
