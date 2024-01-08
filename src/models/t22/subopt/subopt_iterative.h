// Copyright 2023 Eliot Courtney.
#ifndef MODELS_T22_SUBOPT_SUBOPT_ITERATIVE_H_
#define MODELS_T22_SUBOPT_SUBOPT_ITERATIVE_H_

#include <algorithm>
#include <cassert>
#include <compare>
#include <functional>
#include <utility>
#include <vector>

#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "models/t04/mfe/dp.h"
#include "models/t22/energy/model.h"
#include "models/t22/mfe/mfe.h"
#include "models/t22/trace/trace.h"
#include "util/splaymap.h"

namespace mrna::md::t22 {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

class SuboptIterative {
 public:
  SuboptIterative(Primary r, Model::Ptr em, DpState dp, SuboptCfg cfg);

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
  Model::Ptr em_;
  SuboptResult res_;
  DpState dp_;
  SuboptCfg cfg_;

  SplayMap<DpIndex, std::vector<Expansion>> cache_;
  std::vector<Node> q_;
  std::vector<DpIndex> unexpanded_;

  std::pair<int, Energy> RunInternal(
      const SuboptCallback& fn, Energy delta, bool exact_energy, int max);

  const std::vector<Expansion>& GetExpansion(const DpIndex& to_expand) {
    if (!cache_.Find(to_expand)) {
      auto exps = GenerateExpansions(to_expand, cfg_.delta);
      std::sort(exps.begin(), exps.end());
      [[maybe_unused]] auto res = cache_.Insert(to_expand, std::move(exps));
      assert(res);
    }

    return cache_.Get();
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

#endif  // MODELS_T22_SUBOPT_SUBOPT_ITERATIVE_H_
