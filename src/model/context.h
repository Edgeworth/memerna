// Copyright 2016 E.
#ifndef MODEL_CONTEXT_H_
#define MODEL_CONTEXT_H_

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "compute/energy/model.h"
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
  Context(energy::EnergyModel em) : em_(std::move(em)), cfg_() {}
  Context(energy::EnergyModel em, ModelCfg cfg) : em_(std::move(em)), cfg_(std::move(cfg)) {}

  Context() = delete;
  Context(const Context& o) = delete;
  Context& operator=(const Context&) = delete;
  Context(Context&& o) = delete;
  Context& operator=(Context&&) = delete;

  FoldResult Fold(Primary r);
  std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      Primary r, bool sorted, Energy subopt_delta = -1, int subopt_num = -1);
  int Suboptimal(Primary r, subopt::SuboptCallback fn, bool sorted, Energy subopt_delta = -1,
      int subopt_num = -1);
  partition::PartitionResult Partition(Primary r);

 private:
  energy::EnergyModel em_;
  ModelCfg cfg_;

  DpArray ComputeTables(const Primary& r);
};

}  // namespace mrna

#endif  // MODEL_CONTEXT_H_
