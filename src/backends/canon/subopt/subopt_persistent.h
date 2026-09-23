// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_SUBOPT_SUBOPT_PERSISTENT_H_
#define BACKENDS_CANON_SUBOPT_SUBOPT_PERSISTENT_H_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <queue>
#include <string>
#include <tuple>
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
#include "model/base.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::canon {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

template <erg::EnergyCfg CFG, class PseudofreeCfgT, bool UseLru>
class SuboptPersistent {
 public:
  static bool IsSupported(
      const PseudofreeCfgT& pf, const SuboptCfg& subopt_cfg, std::string* reason = nullptr);

  SuboptPersistent(Primary r, Model::Ptr m, DpState dp, PseudofreeCfgT pf, SuboptCfg subopt_cfg);

  using SuboptIndex = uint32_t;
  static constexpr SuboptIndex INVALID_SUBOPT_INDEX = INVALID_UINT32;

  int64_t Run(const SuboptCallback& fn);

 private:
  struct Node {
    SuboptIndex parent_idx{INVALID_SUBOPT_INDEX};
    ExpansionIndex parent_expand_idx{INVALID_EXPANSION_INDEX};
    SuboptIndex unexpanded_idx{INVALID_SUBOPT_INDEX};
    ExpansionIndex unexpanded_expand_idx{INVALID_EXPANSION_INDEX};
    ExpansionIndex expand_idx{0};
    DpIndex to_expand{};
    uint8_t unexpanded_child{0};
  };

  Primary r_;
  Model::Ptr m_;
  PseudofreeCfgT pf_;
  Precomp<CFG, PseudofreeCfgT> pc_;
  DpState dp_;
  SuboptCfg subopt_cfg_;
  SuboptResult res_;

  ExpansionCache<DpIndex, Expansion, UseLru> cache_;
  std::vector<Node> q_;
  std::priority_queue<std::tuple<Energy, SuboptIndex>> pq_;

  std::pair<Energy, SuboptIndex> RunInternal();

  void GenerateResult(SuboptIndex idx);

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

#endif  // BACKENDS_CANON_SUBOPT_SUBOPT_PERSISTENT_H_
