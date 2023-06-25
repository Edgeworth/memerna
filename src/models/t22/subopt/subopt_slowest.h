// Copyright 2023 Eliot Courtney.
#ifndef MODELS_T22_SUBOPT_SUBOPT_FASTEST_H_
#define MODELS_T22_SUBOPT_SUBOPT_FASTEST_H_

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

class SuboptSlowest {
 public:
  SuboptSlowest(Primary r, Model::Ptr em, DpState dp, SuboptCfg cfg);

  int Run(const SuboptCallback& fn);

 private:
  struct DfsState {
    int child_idx = {0};
    std::optional<DpIndex> to_expand = {};
    bool should_unexpand = {false};
  };

  Primary r_;
  Model::Ptr em_;
  SuboptResult res_;
  DpState dp_;
  SuboptCfg cfg_;

  SplayMap<DpIndex, std::vector<Expansion>> cache_;
  std::vector<DfsState> q_;
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

  [[nodiscard]] std::vector<Expansion> ComputeExt(int st, int en, int a);

  [[nodiscard]] std::vector<Expansion> ComputePairedOrNoStack(int st, int en, bool is_nostack);

  [[nodiscard]] std::vector<Expansion> ComputeUnpaired(int st, int en, int a);

  [[nodiscard]] std::vector<Expansion> ComputePenultimate(int st, int en, int length);
};

}  // namespace mrna::md::t22

#endif  // MODELS_T22_SUBOPT_SUBOPT_FASTEST_H_
