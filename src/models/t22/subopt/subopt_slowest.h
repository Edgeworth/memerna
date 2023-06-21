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
#include "util/splaymap.h"

namespace mrna::md::t22 {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

struct Expand {
  // Extra energy of this expansion compared to the best choice.
  Energy energy = {ZERO_E};

  // st is -1 if this does not exist
  DpIndex to_expand = {};
  DpIndex unexpanded = {};
  IndexCtd ctd0 = {};
  IndexCtd ctd1 = {};

  bool operator<(const Expand& o) const { return energy < o.energy; }
};

class SuboptSlowest {
 public:
  SuboptSlowest(Primary r, Model::Ptr em, DpState dp, SuboptCfg cfg);

  int Run(const SuboptCallback& fn);

 private:
  struct DfsState {
    int idx{};
    DpIndex expand;
    bool should_unexpand{};
  };

  Primary r_;
  Model::Ptr em_;
  SuboptResult res_;
  DpState dp_;
  SuboptCfg cfg_;

  SplayMap<DpIndex, std::vector<Expand>> cache_;
  std::vector<DfsState> q_;
  std::vector<DpIndex> unexpanded_;

  std::pair<int, Energy> RunInternal(
      const SuboptCallback& fn, Energy delta, bool exact_energy, int max);

  const std::vector<Expand>& GetExpansion(const DpIndex& to_expand) {
    if (!cache_.Find(to_expand)) {
      auto exps = GenerateExpansions(to_expand, cfg_.delta);
      std::sort(exps.begin(), exps.end());
      [[maybe_unused]] auto res = cache_.Insert(to_expand, std::move(exps));
      assert(res);
    }
    return cache_.Get();
  }

  [[nodiscard]] std::vector<Expand> GenerateExpansions(
      const DpIndex& to_expand, Energy delta) const;
};

}  // namespace mrna::md::t22

#endif  // MODELS_T22_SUBOPT_SUBOPT_FASTEST_H_
