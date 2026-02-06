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
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna {

struct FoldResult {
  mfe::MfeResult mfe;
  trace::TraceResult tb;
};

struct MfeBackend {
  const BackendModelPtr& m;
  MfeAlg alg;
  MfeFn mfe_fn;
  MfeExteriorFn mfe_exterior_fn;
  TraceFn trace_fn;
};

struct SuboptBackend {
  const BackendModelPtr& m;
  SuboptAlg alg;
  MfeAlg mfe_alg;
  SuboptFn subopt_fn;
  MfeFn mfe_fn;
  MfeExteriorFn mfe_exterior_fn;
};

struct PfnBackend {
  const BackendModelPtr& m;
  PfnAlg alg;
  PfnFn pfn_fn;
};

class Ctx {
 public:
  explicit Ctx(std::optional<BackendKind> backend, BackendCfg cfg)
      : cfg_(std::move(cfg)), backend_(backend), backends_() {}
  explicit Ctx(BackendModelPtr m, BackendCfg cfg)
      : cfg_(std::move(cfg)), backend_(GetBackendKind(m)), backends_() {
    backends_[static_cast<int>(*backend_)] = std::move(m);
  }
  ~Ctx() = default;

  Ctx(Ctx&& o) noexcept;
  Ctx& operator=(Ctx&&) noexcept;

  Ctx(const Ctx& o) = delete;
  Ctx& operator=(const Ctx&) = delete;

  [[nodiscard]] erg::EnergyResult Efn(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf = {}, const Ctds* given_ctd = nullptr,
      bool build_structure = false) const;

  [[nodiscard]] FoldResult Fold(const Primary& r, std::optional<MfeAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) const;

  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptIntoVector(const Primary& r,
      std::optional<MfeAlg> mfe_alg, std::optional<SuboptAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, subopt::SuboptCfg subopt_cfg) const;

  [[nodiscard]] int Subopt(const Primary& r, std::optional<MfeAlg> mfe_alg,
      std::optional<SuboptAlg> alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
      const subopt::SuboptCallback& fn, subopt::SuboptCfg subopt_cfg) const;

  [[nodiscard]] pfn::PfnResult Pfn(const Primary& r, std::optional<PfnAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf) const;

  [[nodiscard]] MfeBackend BackendForFold(
      std::optional<MfeAlg> alg, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf) const;
  [[nodiscard]] SuboptBackend BackendForSubopt(std::optional<SuboptAlg> alg,
      std::optional<MfeAlg> mfe_alg, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf,
      const subopt::SuboptCfg& subopt_cfg) const;
  [[nodiscard]] PfnBackend BackendForPfn(
      std::optional<PfnAlg> alg, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf) const;

  static Ctx FromArgParse(const ArgParse& args);

 private:
  BackendCfg cfg_;
  std::optional<BackendKind> backend_;
  mutable std::array<std::optional<BackendModelPtr>, EnumCount<BackendKind>()> backends_;
  mutable std::array<std::once_flag, EnumCount<BackendKind>()> backend_once_;

  [[nodiscard]] const BackendModelPtr& EnsureBackend(BackendKind kind) const;
};

void RegisterOpts(ArgParse* args);

}  // namespace mrna

#endif  // API_CTX_CTX_H_
