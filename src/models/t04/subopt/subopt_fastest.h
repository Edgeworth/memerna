// Copyright 2016 Eliot Courtney.
#ifndef MODELS_T04_SUBOPT_SUBOPT_FASTEST_H_
#define MODELS_T04_SUBOPT_SUBOPT_FASTEST_H_

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
#include "models/t04/energy/model.h"
#include "models/t04/energy/precomp.h"
#include "models/t04/mfe/dp.h"
#include "util/splaymap.h"

namespace mrna::md::t04 {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

struct Expansion {
  // Extra energy of this expansion compared to the best choice.
  Energy delta = {ZERO_E};

  // st is -1 if this does not exist
  DpIndex to_expand = {};
  DpIndex unexpanded = {};
  IndexCtd ctd0 = {};
  IndexCtd ctd1 = {};

  bool operator<(const Expansion& o) const { return delta < o.delta; }
};

class SuboptFastest {
 public:
  SuboptFastest(Primary r, Model::Ptr em, DpState dp, SuboptCfg cfg);

  int Run(const SuboptCallback& fn);

 private:
  struct DfsState {
    // Index of the child expansion of |expand| we should process.
    int child_idx = {0};
    // Index whose child expansions we are processing.
    DpIndex to_expand = {};
    // Stores whether this node's |to_expand| was from |unexpanded_| and needs
    // to be replaced when going back up the DFS stack.
    bool should_unexpand = {false};
  };

  Primary r_;
  Model::Ptr em_;
  Precomp pc_;
  DpState dp_;
  SuboptCfg cfg_;

  SplayMap<DpIndex, std::vector<Expansion>> cache_;
  std::vector<DfsState> q_;

  // Incremental state. Holds the current partial structure.
  SuboptResult res_;
  // Incremental state - holds unexpanded Indexes for the current partial structure.
  std::vector<DpIndex> unexpanded_;

  std::pair<int, Energy> RunInternal(
      const SuboptCallback& fn, Energy delta, bool exact_energy, int max);

  const std::vector<Expansion>& GetExpansion(const DpIndex& to_expand) {
    // We request the expansions of an index multiple times when we find the
    // next sibling of a node after coming back up during the DFS.
    if (!cache_.Find(to_expand)) {
      // Need to generate the full way to delta so we can properly set |next_seen|.
      auto exps = GenerateExpansions(to_expand, cfg_.delta);
      std::sort(exps.begin(), exps.end());
      [[maybe_unused]] auto res = cache_.Insert(to_expand, std::move(exps));
      assert(res);
    }
    return cache_.Get();
  }

  // Generates expansions for the given index, given that the extra energy over the best choice
  // can't be more than |delta|.
  [[nodiscard]] std::vector<Expansion> GenerateExpansions(
      const DpIndex& to_expand, Energy delta) const;
};

}  // namespace mrna::md::t04

#endif  // MODELS_T04_SUBOPT_SUBOPT_FASTEST_H_
