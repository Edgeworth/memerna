// Copyright 2016 Eliot Courtney.
#ifndef CTX_CTX_H_
#define CTX_CTX_H_

#include <memory>
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
  explicit Ctx(energy::EnergyModelPtr em) : em_(std::move(em)), cfg_() {}
  ~Ctx() = default;
  Ctx(energy::EnergyModelPtr em, CtxCfg cfg) : em_(std::move(em)), cfg_(cfg) {}

  Ctx(Ctx&& o) = default;
  Ctx& operator=(Ctx&&) = default;

  Ctx(const Ctx& o) = delete;
  Ctx& operator=(const Ctx&) = delete;

  energy::EnergyResult Efn(const Primary& r, const Secondary& s, const Ctds* given_ctd = nullptr,
      bool build_structure = false) const;
  [[nodiscard]] FoldResult Fold(const Primary& r) const;
  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      const Primary& r, subopt::SuboptCfg cfg) const;
  [[nodiscard]] int Suboptimal(
      const Primary& r, const subopt::SuboptCallback& fn, subopt::SuboptCfg cfg) const;
  [[nodiscard]] part::PartResult Partition(const Primary& r) const;

  [[nodiscard]] const energy::EnergyModel& em() const { return *em_; }

  static Ctx FromArgParse(const ArgParse& args);

 private:
  energy::EnergyModelPtr em_;
  CtxCfg cfg_;

  [[nodiscard]] DpArray ComputeTables(const Primary& r) const;
};

}  // namespace mrna::ctx

#endif  // CTX_CTX_H_
