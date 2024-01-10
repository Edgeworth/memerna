// Copyright 2023 Eliot Courtney.

#ifndef MODELS_T22_SUBOPT_SUBOPT_PERSISTENT_H_
#define MODELS_T22_SUBOPT_SUBOPT_PERSISTENT_H_

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
#include "models/t22/energy/model.h"
#include "models/t22/trace/trace.h"
#include "util/splaymap.h"

namespace mrna::md::t22 {

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
    std::optional<DpIndex> to_expand{};
  };

  Primary r_;
  Model::Ptr em_;
  DpState dp_;
  SuboptCfg cfg_;
  SuboptResult res_;

  std::vector<std::vector<Expansion>> cache_;
  std::vector<Node> q_;
  std::priority_queue<std::pair<Energy, int>> pq_;

  std::pair<Energy, int> RunInternal();

  // Computes the suboptimal folding for the given subpath and puts it into `res_`.
  void GenerateResult(int idx);

  const std::vector<Expansion>& GetExpansion(const DpIndex& to_expand) {
    auto idx = LinearIndex(to_expand, r_.size());
    if (cache_[idx].empty()) {
      auto exps = GenerateExpansions(to_expand, cfg_.delta);
      std::sort(exps.begin(), exps.end());
      assert(!exps.empty());
      cache_[idx] = std::move(exps);
    }
    return cache_[idx];
  }

  [[nodiscard]] std::vector<Expansion> GenerateExpansions(
      const DpIndex& to_expand, Energy delta) const;

  [[nodiscard]] std::vector<Expansion> ExtExpansions(int st, int a, Energy delta) const;

  [[nodiscard]] std::vector<Expansion> PairedOrNoStackExpansions(
      int st, int en, bool is_nostack, Energy delta) const;

  [[nodiscard]] std::vector<Expansion> UnpairedExpansions(
      int st, int en, int a, Energy delta) const;

  [[nodiscard]] std::vector<Expansion> PenultimateExpansions(
      int st, int en, int length, Energy delta) const;
};

}  // namespace mrna::md::t22

#endif  // MODELS_T22_SUBOPT_SUBOPT_PERSISTENT_H_
