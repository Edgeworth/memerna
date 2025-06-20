// Copyright 2016 E.
#ifndef BACKENDS_BASE_SUBOPT_SUBOPT_DEBUG_H_
#define BACKENDS_BASE_SUBOPT_SUBOPT_DEBUG_H_

#include <set>
#include <vector>

#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

using mrna::subopt::SuboptCallback;
using mrna::subopt::SuboptCfg;
using mrna::subopt::SuboptResult;

class SuboptDebug {
 public:
  SuboptDebug(Primary r, Model::Ptr m, DpState dp, SuboptCfg cfg);

  int Run(const SuboptCallback& fn);

 private:
  struct Node {
    // State should be fully defined by `not_yet_expanded` and `ctd` which denote
    // what it has done so far, and what it can do from now.
    std::vector<DpIndex> not_yet_expanded;
    SuboptResult res;  // Stores the minimum energy this state could have.

    [[nodiscard]] Node copy() const {
      return Node{.not_yet_expanded = not_yet_expanded, .res = SuboptResult(res)};
    }

    bool operator<(const Node& o) const { return res.energy < o.res.energy; }
  };

  Primary r_;
  Model::Ptr m_;
  DpState dp_;
  SuboptCfg cfg_;

  // This node is where we build intermediate results to be pushed onto the queue.
  Node curnode_;
  std::multiset<Node> finished_;
  std::multiset<Node> q_;

  void PruneInsert(const Node& node, std::multiset<Node>* prune) {
    if (node.res.energy <= dp_.ext[0][EXT] + cfg_.delta) {
      if (static_cast<int>(prune->size()) >= cfg_.strucs &&
          (--prune->end())->res.energy > node.res.energy)
        prune->erase(--prune->end());
      if (static_cast<int>(prune->size()) < cfg_.strucs) prune->insert(node.copy());
    }
  }

  // Creates and inserts a new node with energy `energy` that doesn't
  // need to expand any more ranges than it currently has.
  void Expand(Energy energy) {
    curnode_.res.energy = energy;
    PruneInsert(curnode_, &q_);
  }

  // Creates and inserts a new node with energy `energy` that needs to expand the given ranges.
  void Expand(Energy energy, DpIndex nye) {
    curnode_.not_yet_expanded.push_back(nye);
    curnode_.res.energy = energy;
    PruneInsert(curnode_, &q_);
    curnode_.not_yet_expanded.pop_back();
  }

  void Expand(Energy energy, DpIndex nye, IndexCtd ctd_idx) {
    curnode_.res.tb.ctd[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye);
    curnode_.res.tb.ctd[ctd_idx.idx] = CTD_NA;
  }

  // Creates and inserts a new node with energy `energy` that needs to expand the two given ranges.
  void Expand(Energy energy, DpIndex nye0, DpIndex nye1) {
    curnode_.not_yet_expanded.push_back(nye0);
    curnode_.not_yet_expanded.push_back(nye1);
    curnode_.res.energy = energy;
    PruneInsert(curnode_, &q_);
    curnode_.not_yet_expanded.pop_back();
    curnode_.not_yet_expanded.pop_back();
  }

  void Expand(Energy energy, DpIndex nye0, DpIndex nye1, IndexCtd ctd_idx) {
    curnode_.res.tb.ctd[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye0, nye1);
    curnode_.res.tb.ctd[ctd_idx.idx] = CTD_NA;
  }

  void Expand(Energy energy, DpIndex nye0, DpIndex nye1, IndexCtd ctd_idx0, IndexCtd ctd_idx1) {
    curnode_.res.tb.ctd[ctd_idx0.idx] = ctd_idx0.ctd;
    curnode_.res.tb.ctd[ctd_idx1.idx] = ctd_idx1.ctd;
    Expand(energy, nye0, nye1);
    curnode_.res.tb.ctd[ctd_idx0.idx] = CTD_NA;
    curnode_.res.tb.ctd[ctd_idx1.idx] = CTD_NA;
  }
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_SUBOPT_SUBOPT_DEBUG_H_
