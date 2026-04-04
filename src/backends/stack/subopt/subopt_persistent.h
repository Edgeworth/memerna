// Copyright 2023 Eliot Courtney.

#ifndef BACKENDS_STACK_SUBOPT_SUBOPT_PERSISTENT_H_
#define BACKENDS_STACK_SUBOPT_SUBOPT_PERSISTENT_H_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "backends/common/energy.h"
#include "backends/common/expansion_cache.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/trace/trace.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::stack {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

// Subopt folding based on a persistent data structure algorithm.
template <bool UseLru>
class SuboptPersistent {
 public:
  static bool IsSupported(const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf,
      const SuboptCfg& subopt_cfg, std::string* reason = nullptr);

  SuboptPersistent(Primary r, Model::Ptr m, DpState dp, erg::EnergyCfg cfg, erg::PseudofreeCfg pf,
      SuboptCfg subopt_cfg);

  using SuboptIndex = uint32_t;
  static constexpr SuboptIndex INVALID_SUBOPT_INDEX = INVALID_UINT32;

  int64_t Run(const SuboptCallback& fn);

 private:
  struct Node {
    // Index of the parent Node in the expand tree.
    SuboptIndex parent_idx{INVALID_SUBOPT_INDEX};
    // Index of the expansion this Node used w.r.t. the parent state's `to_expand`.
    // INVALID_EXPANSION_INDEX means root.
    ExpansionIndex parent_expand_idx{INVALID_EXPANSION_INDEX};
    // Index of the next Node who's expansion contains an unexpanded DpIndex we need to process.
    SuboptIndex unexpanded_idx{INVALID_SUBOPT_INDEX};
    // Index of the expansion to use for `unexpanded_idx`'s Node. INVALID_EXPANSION_INDEX means
    // none.
    ExpansionIndex unexpanded_expand_idx{INVALID_EXPANSION_INDEX};
    // Index of the child expansion of `to_expand` we should process. This gets updated in place
    // (saves time and memory), which is why we need to keep `parent_expand_idx` around as well.
    ExpansionIndex expand_idx{0};
    // DpIndex whose child expansions we are processing.
    std::optional<DpIndex> to_expand{};
  };

  Primary r_;
  Model::Ptr m_;
  erg::EnergyCfg cfg_;
  erg::PseudofreeCfg pf_;
  DpState dp_;
  SuboptCfg subopt_cfg_;
  SuboptResult res_;

  ExpansionCache<DpIndex, Expansion, UseLru> cache_;
  std::vector<Node> q_;
  std::priority_queue<std::pair<Energy, SuboptIndex>> pq_;

  std::pair<Energy, SuboptIndex> RunInternal();

  // Computes the suboptimal folding for the given subpath and puts it into `res_`.
  void GenerateResult(SuboptIndex idx);

  const std::vector<Expansion>& GetExpansion(const DpIndex& to_expand) {
    auto key = LinearIndex(to_expand, r_.size());
    if (const auto& val = cache_.Get(key); !val.empty()) return val;

    auto exps = GenerateExpansions(to_expand, subopt_cfg_.delta);
    std::sort(exps.begin(), exps.end());
    assert(!exps.empty());
    return cache_.Insert(key, std::move(exps));
  }

  [[nodiscard]] std::vector<Expansion> GenerateExpansions(
      const DpIndex& to_expand, Energy delta) const;

  [[nodiscard]] std::vector<Expansion> ExtExpansions(int st, base::DpArrayId a, Energy delta) const;

  [[nodiscard]] std::vector<Expansion> PairedOrNoStackExpansions(
      int st, int en, bool is_nostack, Energy delta) const;

  [[nodiscard]] std::vector<Expansion> UnpairedExpansions(
      int st, int en, base::DpArrayId a, Energy delta) const;

  [[nodiscard]] std::vector<Expansion> PenultimateExpansions(
      int st, int en, int length, Energy delta) const;
};

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_SUBOPT_SUBOPT_PERSISTENT_H_
