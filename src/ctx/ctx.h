// Copyright 2016 E.
#ifndef CTX_CTX_H_
#define CTX_CTX_H_

#include <utility>
#include <vector>

#include "compute/dp.h"
#include "compute/energy/energy.h"
#include "compute/energy/model.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/partition.h"
#include "compute/subopt/config.h"
#include "compute/subopt/subopt.h"
#include "compute/traceback/traceback.h"
#include "ctx/config.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::ctx {

struct FoldResult {
  mfe::MfeResult mfe;
  tb::TracebackResult tb;
};

class Ctx {
 public:
  explicit Ctx(energy::EnergyModel em) : em_(std::move(em)), cfg_() {}
  Ctx(energy::EnergyModel em, CtxCfg cfg) : em_(std::move(em)), cfg_(std::move(cfg)) {}

  Ctx(Ctx&& o) = default;
  Ctx& operator=(Ctx&&) = default;

  Ctx(const Ctx& o) = delete;
  Ctx& operator=(const Ctx&) = delete;

  energy::EnergyResult Efn(
      Primary r, Secondary s, const Ctds* given_ctd = nullptr, bool build_structure = false) const;
  FoldResult Fold(Primary r) const;
  std::vector<subopt::SuboptResult> SuboptimalIntoVector(Primary r, subopt::SuboptCfg cfg) const;
  int Suboptimal(Primary r, subopt::SuboptCallback fn, subopt::SuboptCfg cfg) const;
  part::PartResult Partition(Primary r) const;

  const energy::EnergyModel& em() const { return em_; }

  static Ctx FromArgParse(const ArgParse& args);

 private:
  energy::EnergyModel em_;
  CtxCfg cfg_;

  DpArray ComputeTables(const Primary& r) const;
};

}  // namespace mrna::ctx

#endif  // CTX_CTX_H_
