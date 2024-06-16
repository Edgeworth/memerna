// Copyright 2016 Eliot Courtney.
#ifndef API_CTX_CTX_H_
#define API_CTX_CTX_H_

#include <utility>
#include <vector>

#include "api/ctx/ctx_cfg.h"
#include "api/energy/energy.h"
#include "api/mfe.h"
#include "api/pfn.h"
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
  explicit Ctx(BackendModelPtr m) : m_(std::move(m)), cfg_() {}
  ~Ctx() = default;
  Ctx(BackendModelPtr m, CtxCfg cfg) : m_(std::move(m)), cfg_(cfg) {}

  Ctx(Ctx&& o) = default;
  Ctx& operator=(Ctx&&) = default;

  Ctx(const Ctx& o) = delete;
  Ctx& operator=(const Ctx&) = delete;

  erg::EnergyResult Efn(const Primary& r, const Secondary& s, const Ctds* given_ctd = nullptr,
      bool build_structure = false) const;
  [[nodiscard]] FoldResult Fold(const Primary& r, const trace::TraceCfg& cfg) const;
  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptIntoVector(
      const Primary& r, subopt::SuboptCfg cfg) const;
  [[nodiscard]] int Subopt(
      const Primary& r, const subopt::SuboptCallback& fn, subopt::SuboptCfg cfg) const;
  [[nodiscard]] pfn::PfnResult Pfn(const Primary& r) const;

  [[nodiscard]] const BackendModelPtr& m() const { return m_; }

  static Ctx FromArgParse(const ArgParse& args);

 private:
  BackendModelPtr m_;
  CtxCfg cfg_;

  void ComputeMfe(const Primary& r, mfe::DpState& dp) const;
  Energy ComputeMfeExterior(const Primary& r, mfe::DpState& dp) const;
  [[nodiscard]] trace::TraceResult ComputeTraceback(
      const Primary& r, const trace::TraceCfg& cfg, const mfe::DpState& dp) const;
};

}  // namespace mrna

#endif  // API_CTX_CTX_H_
