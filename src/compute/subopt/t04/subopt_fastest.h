// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_SUBOPT_T04_SUBOPT_FASTEST_H_
#define COMPUTE_SUBOPT_T04_SUBOPT_FASTEST_H_

#include <algorithm>
#include <cassert>
#include <compare>
#include <utility>
#include <vector>

#include "compute/energy/t04/model.h"
#include "compute/energy/t04/precomp.h"
#include "compute/mfe/t04/dp.h"
#include "compute/subopt/subopt.h"
#include "compute/subopt/subopt_cfg.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/splaymap.h"

namespace mrna::subopt::t04 {

struct Expand {
  Expand() = delete;
  explicit Expand(Energy energy_) : energy(energy_) {}
  Expand(Energy energy_, const Index& to_expand_) : energy(energy_), to_expand(to_expand_) {}
  Expand(Energy energy_, const Index& to_expand_, const Index& unexpanded_)
      : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_) {}
  Expand(Energy energy_, const Index& to_expand_, const IndexCtd& ctd0_)
      : energy(energy_), to_expand(to_expand_), ctd0(ctd0_) {}
  Expand(Energy energy_, const Index& to_expand_, const Index& unexpanded_, const IndexCtd& ctd0_)
      : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_), ctd0(ctd0_) {}
  Expand(Energy energy_, const Index& to_expand_, const Index& unexpanded_, const IndexCtd& ctd0_,
      const IndexCtd& ctd1_)
      : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_), ctd0(ctd0_), ctd1(ctd1_) {}
  Energy energy;
  Index to_expand;  // st is -1 if this does not exist
  Index unexpanded;
  IndexCtd ctd0;
  IndexCtd ctd1;

  bool operator<(const Expand& o) const { return energy < o.energy; }
};

class SuboptFastest {
 public:
  SuboptFastest(Primary r, erg::t04::Model::Ptr em, DpArray dp, ExtArray ext, SuboptCfg cfg);

  int Run(const SuboptCallback& fn);

 private:
  struct DfsState {
    int idx{};
    Index expand;
    // Stores whether this node's |expand| was from |unexpanded| and needs to be replaced.
    bool should_unexpand{};
  };

  Primary r_;
  erg::t04::Model::Ptr em_;
  erg::t04::Precomp pc_;
  SuboptResult res_;
  DpArray dp_;
  ExtArray ext_;
  SuboptCfg cfg_;

  SplayMap<Index, std::vector<Expand>> cache_;
  std::vector<DfsState> q_;
  std::vector<Index> unexpanded_;

  std::pair<int, Energy> RunInternal(
      const SuboptCallback& fn, Energy delta, bool exact_energy, int max);

  const std::vector<Expand>& GetExpansion(const Index& to_expand) {
    if (!cache_.Find(to_expand)) {
      // Need to generate the full way to delta so we can properly set |next_seen|.
      auto exps = GenerateExpansions(to_expand, cfg_.delta);
      std::sort(exps.begin(), exps.end());
      [[maybe_unused]] auto res = cache_.Insert(to_expand, std::move(exps));
      assert(res);
    }
    return cache_.Get();
  }

  [[nodiscard]] std::vector<Expand> GenerateExpansions(const Index& to_expand, Energy delta) const;
};

}  // namespace mrna::subopt::t04

#endif  // COMPUTE_SUBOPT_T04_SUBOPT_FASTEST_H_
