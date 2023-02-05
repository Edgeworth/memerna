// Copyright 2016 E.
#ifndef COMPUTE_SUBOPT_T04_SUBOPT_SLOWEST_H_
#define COMPUTE_SUBOPT_T04_SUBOPT_SLOWEST_H_

#include <compare>
#include <set>
#include <vector>

#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "compute/traceback/traceback.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "models/t04/energy/model.h"
#include "models/t04/mfe/dp.h"

namespace mrna::md::t04::subopt {

class SuboptSlowest {
 public:
  SuboptSlowest(Primary r, erg::t04::Model::Ptr em, mfe::t04::DpState dp, SuboptCfg cfg);

  int Run(const SuboptCallback& fn);

 private:
  struct Node {
    // State should be fully defined by |not_yet_expanded|, |history|, and |ctd| which denote
    // what it has done so far, and what it can do from now.
    std::vector<Index> not_yet_expanded;
    std::vector<Index> history;
    SuboptResult res;  // Stores the minimum energy this state could have.

    [[nodiscard]] Node copy() const { return Node{not_yet_expanded, history, SuboptResult(res)}; }

    bool operator<(const Node& o) const { return res.energy < o.res.energy; }
  };

  Primary r_;
  erg::t04::Model::Ptr em_;
  DpArray dp_;
  ExtArray ext_;
  SuboptCfg cfg_;

  // This node is where we build intermediate results to be pushed onto the queue.
  Node curnode_;
  std::multiset<Node> finished_;
  std::multiset<Node> q_;

  void PruneInsert(const Node& node, std::multiset<Node>* prune) {
    if (node.res.energy <= ext_[0][EXT] + cfg_.delta) {
      if (static_cast<int>(prune->size()) >= cfg_.strucs &&
          (--prune->end())->res.energy > node.res.energy)
        prune->erase(--prune->end());
      if (static_cast<int>(prune->size()) < cfg_.strucs) prune->insert(node.copy());
    }
  }

  // Creates and inserts a new node with energy |energy| that doesn't
  // need to expand any more ranges than it currently has.
  void Expand(Energy energy) {
    curnode_.res.energy = energy;
    PruneInsert(curnode_, &q_);
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the given ranges.
  void Expand(Energy energy, Index nye) {
    curnode_.not_yet_expanded.push_back(nye);
    curnode_.res.energy = energy;
    PruneInsert(curnode_, &q_);
    curnode_.not_yet_expanded.pop_back();
  }

  void Expand(Energy energy, Index nye, IndexCtd ctd_idx) {
    curnode_.res.tb.ctd[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye);
    curnode_.res.tb.ctd[ctd_idx.idx] = CTD_NA;
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the two given ranges.
  void Expand(Energy energy, Index nye0, Index nye1) {
    curnode_.not_yet_expanded.push_back(nye0);
    curnode_.not_yet_expanded.push_back(nye1);
    curnode_.res.energy = energy;
    PruneInsert(curnode_, &q_);
    curnode_.not_yet_expanded.pop_back();
    curnode_.not_yet_expanded.pop_back();
  }

  void Expand(Energy energy, Index nye0, Index nye1, IndexCtd ctd_idx) {
    curnode_.res.tb.ctd[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye0, nye1);
    curnode_.res.tb.ctd[ctd_idx.idx] = CTD_NA;
  }

  void Expand(Energy energy, Index nye0, Index nye1, IndexCtd ctd_idx0, IndexCtd ctd_idx1) {
    curnode_.res.tb.ctd[ctd_idx0.idx] = ctd_idx0.ctd;
    curnode_.res.tb.ctd[ctd_idx1.idx] = ctd_idx1.ctd;
    Expand(energy, nye0, nye1);
    curnode_.res.tb.ctd[ctd_idx0.idx] = CTD_NA;
    curnode_.res.tb.ctd[ctd_idx1.idx] = CTD_NA;
  }
};

}  // namespace mrna::md::t04::subopt

#endif  // COMPUTE_SUBOPT_T04_SUBOPT_SLOWEST_H_