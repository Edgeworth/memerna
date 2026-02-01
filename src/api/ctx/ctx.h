// Copyright 2016 Eliot Courtney.
#ifndef API_CTX_CTX_H_
#define API_CTX_CTX_H_

#include <array>
#include <mutex>
#include <optional>
#include <utility>
#include <vector>

#include "api/ctx/algorithm.h"
#include "api/ctx/backend.h"
#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
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
  explicit Ctx(BackendCfg cfg) : cfg_(std::move(cfg)), backends_() {}
  explicit Ctx(BackendModelPtr m) : cfg_(std::nullopt), backends_() { backends_[0] = std::move(m); }
  ~Ctx() = default;

  Ctx(Ctx&& o) noexcept;
  Ctx& operator=(Ctx&&) noexcept;

  Ctx(const Ctx& o) = delete;
  Ctx& operator=(const Ctx&) = delete;

  erg::EnergyResult Efn(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf = {}, const Ctds* given_ctd = nullptr,
      bool build_structure = false) const;

  [[nodiscard]] FoldResult Fold(const Primary& r, MfeAlg alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) const;

  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptIntoVector(const Primary& r, MfeAlg mfe_alg,
      SuboptAlg alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
      subopt::SuboptCfg subopt_cfg) const;

  [[nodiscard]] int Subopt(const Primary& r, MfeAlg mfe_alg, SuboptAlg alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
      subopt::SuboptCfg subopt_cfg) const;

  [[nodiscard]] pfn::PfnResult Pfn(
      const Primary& r, PfnAlg alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf) const;

  // Access the backend model (loads if necessary)
  [[nodiscard]] const BackendModelPtr& m() const { return EnsureBackend(); }

  static Ctx FromArgParse(const ArgParse& args);

 private:
  std::optional<BackendCfg> cfg_;
  mutable std::array<std::optional<BackendModelPtr>, EnumCount<BackendKind>()> backends_;
  mutable std::array<std::once_flag, EnumCount<BackendKind>()> backend_once_;

  [[nodiscard]] const BackendModelPtr& EnsureBackend() const;

  void ComputeMfe(const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, MfeAlg alg,
      erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf) const;
  Energy ComputeMfeExterior(const BackendModelPtr& m, const Primary& r, mfe::DpState& dp,
      erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf) const;
  [[nodiscard]] trace::TraceResult ComputeTraceback(const BackendModelPtr& m, const Primary& r,
      const mfe::DpState& dp, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
      const trace::TraceCfg& trace_cfg) const;
};

void RegisterOpts(ArgParse* args);

}  // namespace mrna

#endif  // API_CTX_CTX_H_
