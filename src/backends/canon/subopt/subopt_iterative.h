// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_SUBOPT_SUBOPT_ITERATIVE_H_
#define BACKENDS_CANON_SUBOPT_SUBOPT_ITERATIVE_H_

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "backends/canon/dp.h"
#include "backends/canon/energy/model.h"
#include "backends/canon/energy/precomp.h"
#include "backends/common/energy.h"
#include "backends/common/expansion_cache.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::canon {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

template <erg::EnergyCfg CFG, class PseudofreeCfgT, bool UseLru>
class SuboptIterative {
 public:
  static bool IsSupported(
      const PseudofreeCfgT& pf, const SuboptCfg& subopt_cfg, std::string* reason = nullptr);

  SuboptIterative(Primary r, Model::Ptr m, DpState dp, PseudofreeCfgT pf, SuboptCfg subopt_cfg);

  int64_t Run(const SuboptCallback& fn);

 private:
  struct Node {
    ExpansionIndex expand_idx = {0};
    DpIndex to_expand{};
    bool should_unexpand = {false};
  };

  Primary r_;
  Model::Ptr m_;
  PseudofreeCfgT pf_;
  Precomp<CFG, PseudofreeCfgT> pc_;
  DpState dp_;
  SuboptCfg subopt_cfg_;

  ExpansionCache<DpIndex, Expansion, UseLru> cache_;
  std::vector<Node> q_;

  SuboptResult res_;
  std::vector<DpIndex> unexpanded_;
  int iters_until_check_ = 1;

  struct RunInternalResult {
    int64_t count = 0;
    Energy next_seen = MAX_E;
    bool keep_going = true;
  };

  template <bool CheckTime>
  RunInternalResult RunInternal(const SuboptCallback& fn, Energy delta, bool exact_energy,
      int64_t max, const std::chrono::steady_clock::time_point& start_time);

  const std::vector<Expansion>& GetExpansion(const DpIndex& to_expand) {
    auto key = to_expand.LinearIndex(r_.size());
    if (const auto& val = cache_.Get(key); !val.empty()) return val;

    auto exps = GenerateExpansions(to_expand, subopt_cfg_.delta);
    std::sort(exps.begin(), exps.end());
    assert(!exps.empty());
    return cache_.Insert(key, std::move(exps));
  }

  [[nodiscard]] std::vector<Expansion> GenerateExpansions(
      const DpIndex& to_expand, Energy delta) const;
};

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_SUBOPT_SUBOPT_ITERATIVE_H_
