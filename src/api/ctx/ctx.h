// Copyright 2016 Eliot Courtney.
#ifndef API_CTX_CTX_H_
#define API_CTX_CTX_H_

#include <utility>
#include <vector>

#include "api/ctx/ctx_cfg.h"
#include "api/energy/energy.h"
#include "api/energy/model.h"
#include "api/mfe.h"
#include "api/part.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "api/trace/trace_cfg.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna {

struct FoldResult {
  mfe::MfeResult mfe;
  trace::TraceResult tb;
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
  [[nodiscard]] FoldResult Fold(const Primary& r, const trace::TraceCfg& cfg) const;
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

  void ComputeMfe(const Primary& r, mfe::DpState& dp) const;
  Energy ComputeMfeExterior(const Primary& r, mfe::DpState& dp) const;
  [[nodiscard]] trace::TraceResult ComputeTraceback(
      const Primary& r, const trace::TraceCfg& cfg, const mfe::DpState& dp) const;
};

}  // namespace mrna

#endif  // API_CTX_CTX_H_
