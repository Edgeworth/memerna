// Copyright 2016 E.
#ifndef COMPUTE_SUBOPT_SUBOPT1_H_
#define COMPUTE_SUBOPT_SUBOPT1_H_

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "compute/energy/precomp.h"
#include "compute/subopt/subopt.h"
#include "model/model.h"
#include "model/primary.h"
#include "util/splaymap.h"

namespace mrna::subopt {

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

class Suboptimal1 {
 public:
  Suboptimal1(Primary r, energy::EnergyModel em, DpArray dp, ExtArray ext, Energy delta, int num);

  int Run(SuboptCallback fn, bool sorted);

 private:
  struct DfsState {
    int idx;
    Index expand;
    // Stores whether this node's |expand| was from |unexpanded| and needs to be replaced.
    bool should_unexpand;
  };

  Primary r_;
  energy::EnergyModel em_;
  energy::Precomp pc_;
  SuboptResult res_;
  DpArray dp_;
  ExtArray ext_;

  const Energy delta_;
  const int max_structures_;
  SplayMap<Index, std::vector<Expand>> cache_;
  std::vector<DfsState> q_;
  std::vector<Index> unexpanded_;

  std::pair<int, int> RunInternal(
      SuboptCallback fn, Energy cur_delta, bool exact_energy, int structure_limit);

  const std::vector<Expand>& GetExpansion(const Index& to_expand) {
    if (!cache_.Find(to_expand)) {
      // Need to generate the full way to delta so we can properly set |next_seen|.
      auto exps = GenerateExpansions(to_expand, delta_);
      std::sort(exps.begin(), exps.end());
      [[maybe_unused]] auto res = cache_.Insert(to_expand, std::move(exps));
      assert(res);
    }
    return cache_.Get();
  }

  std::vector<Expand> GenerateExpansions(const Index& to_expand, Energy delta) const;
};

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_SUBOPT1_H_
