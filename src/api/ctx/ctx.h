// Copyright 2016 E.
#ifndef CTX_CTX_H_
#define CTX_CTX_H_

#include <utility>
#include <vector>

#include "api/energy/model.h"
#include "api/part.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "compute/mfe/mfe.h"
#include "compute/traceback/traceback.h"
#include "ctx/ctx_cfg.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/t04/mfe/dp.h"
#include "util/argparse.h"

namespace mrna {

struct FoldResult {
  mfe::MfeResult mfe;
  tb::TraceResult tb;
};

class Ctx {
 public:
  explicit Ctx(erg::EnergyModelPtr em) : em_(std::move(em)), cfg_() {}
  ~Ctx() = default;
  Ctx(erg::EnergyModelPtr em, CtxCfg cfg) : em_(std::move(em)), cfg_(cfg) {}

  Ctx(Ctx&& o) = default;
  Ctx& operator=(Ctx&&) = default;

  Ctx(const Ctx& o) = delete;
  Ctx& operator=(const Ctx&) = delete;

  erg::EnergyResult Efn(const Primary& r, const Secondary& s, const Ctds* given_ctd = nullptr,
      bool build_structure = false) const;
  [[nodiscard]] FoldResult Fold(const Primary& r) const;
  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      const Primary& r, subopt::SuboptCfg cfg) const;
  [[nodiscard]] int Suboptimal(
      const Primary& r, const subopt::SuboptCallback& fn, subopt::SuboptCfg cfg) const;
  [[nodiscard]] part::PartResult Partition(const Primary& r) const;

  [[nodiscard]] const erg::EnergyModelPtr& em() const { return em_; }

  static Ctx FromArgParse(const ArgParse& args);

 private:
  erg::EnergyModelPtr em_;
  CtxCfg cfg_;

  [[nodiscard]] DpArray ComputeMfe(const Primary& r) const;
  [[nodiscard]] ExtArray ComputeMfeExterior(const Primary& r, const DpArray& dp) const;
  [[nodiscard]] tb::TraceResult ComputeTraceback(
      const Primary& r, const DpArray& dp, const ExtArray& ext) const;
};

}  // namespace mrna

#endif  // CTX_CTX_H_