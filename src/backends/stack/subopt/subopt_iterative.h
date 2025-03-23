// Copyright 2023 Eliot Courtney.
#ifndef BACKENDS_STACK_SUBOPT_SUBOPT_ITERATIVE_H_
#define BACKENDS_STACK_SUBOPT_SUBOPT_ITERATIVE_H_

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "backends/common/expansion_cache.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/mfe.h"
#include "backends/stack/trace/trace.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::stack {

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
    int expand_idx = {0};
    // We could use std::monostate in DpIndex to avoid optional here, but it doesn't change the size
    // of `Node` and it causes a perf regression by adding more variants which worsens the
    // `SplayMap` cache lookup.
    std::optional<DpIndex> to_expand = {};
    bool should_unexpand = {false};
  };

  Primary r_;
  Model::Ptr m_;
  SuboptResult res_;
  DpState dp_;
  SuboptCfg cfg_;

  ExpansionCache<DpIndex, Expansion, UseLru> cache_;
  std::vector<Node> q_;
  std::vector<DpIndex> unexpanded_;

  std::pair<int, Energy> RunInternal(
      const SuboptCallback& fn, Energy delta, bool exact_energy, int max);

  const std::vector<Expansion>& GetExpansion(const DpIndex& to_expand) {
    auto key = LinearIndex(to_expand, r_.size());
    if (const auto& val = cache_.Get(key); !val.empty()) return val;

    auto exps = GenerateExpansions(to_expand, cfg_.delta);
    std::sort(exps.begin(), exps.end());
    assert(!exps.empty());
    return cache_.Insert(key, std::move(exps));
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

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_SUBOPT_SUBOPT_ITERATIVE_H_
