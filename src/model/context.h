// Copyright 2016 Eliot Courtney.
#ifndef MODEL_CONTEXT_H_
#define MODEL_CONTEXT_H_

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "compute/energy/energy.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "compute/traceback/traceback.h"
#include "model/config.h"
#include "model/model.h"
#include "util/argparse.h"

namespace mrna {

struct FoldResult {
  mfe::MfeResult mfe;
  tb::TracebackResult tb;
};

class Context {
 public:
  Context(Primary r, energy::EnergyModel em) : r_(std::move(r)), em_(std::move(em)), cfg_() {
    verify(r_.size() > 0u, "cannot process zero length RNA");
  }
  Context(Primary r, energy::EnergyModel em, ModelCfg cfg)
      : r_(std::move(r)), em_(std::move(em)), cfg_(std::move(cfg)) {
    verify(r_.size() > 0u, "cannot process zero length RNA");
  }
  Context(const Context& o) = default;

  Context() = delete;
  Context& operator=(const Context&) = delete;
  Context(Context&& o) = delete;
  Context& operator=(Context&&) = delete;

  FoldResult Fold();
  std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      bool sorted, Energy subopt_delta = -1, int subopt_num = -1);
  int Suboptimal(
      subopt::SuboptCallback fn, bool sorted, Energy subopt_delta = -1, int subopt_num = -1);
  partition::PartitionResult Partition();

 private:
  Primary r_;
  energy::EnergyModel em_;
  ModelCfg cfg_;

  DpArray ComputeTables();
};

}  // namespace mrna

#endif  // MODEL_CONTEXT_H_
