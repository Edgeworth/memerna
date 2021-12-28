// Copyright 2016 E.
#ifndef MODEL_CONTEXT_H_
#define MODEL_CONTEXT_H_

#include <map>
#include <string>
#include <vector>

#include "compute/energy/energy.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "model/config.h"
#include "model/structure.h"
#include "util/argparse.h"

namespace mrna {

class Context {
 public:
  Context(const Primary& r, const energy::EnergyModelPtr em) : r(r), em(em), cfg() {
    verify(r.size() > 0u, "cannot process zero length RNA");
  }
  Context(const Primary& r, const energy::EnergyModelPtr em, ModelCfg cfg)
      : r(r), em(em), cfg(cfg) {
    verify(r.size() > 0u, "cannot process zero length RNA");
  }
  Context(const Context& o) = default;

  Context() = delete;
  Context& operator=(const Context&) = delete;
  Context(Context&& o) = delete;
  Context& operator=(Context&&) = delete;

  Computed Fold();
  std::vector<Computed> SuboptimalIntoVector(
      bool sorted, Energy subopt_delta = -1, int subopt_num = -1);
  int Suboptimal(
      subopt::SuboptimalCallback fn, bool sorted, Energy subopt_delta = -1, int subopt_num = -1);
  partition::Partition Partition();

 private:
  const Primary r;
  const energy::EnergyModelPtr em;
  const ModelCfg cfg;

  void ComputeTables();
};

}  // namespace mrna

#endif  // MODEL_CONTEXT_H_
