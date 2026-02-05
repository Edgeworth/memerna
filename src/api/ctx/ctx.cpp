// Copyright 2016 Eliot Courtney.
#include "api/ctx/ctx.h"

#include <cassert>
#include <new>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "api/ctx/backend.h"
#include "api/mfe.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace_cfg.h"
#include "backends/brute/alg.h"
#include "backends/common/base/dp.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/dp.h"
#include "model/primary.h"
#include "util/error.h"
#include "util/util.h"

namespace mrna {

namespace {

mfe::DpState CreateDpState(const BackendModelPtr& m) {
  auto vis = overloaded{
      [&](const md::base::Model::Ptr&) -> mfe::DpState { return md::base::DpState{}; },
      [&](const md::base::opt::Model::Ptr&) -> mfe::DpState { return md::base::DpState{}; },
      [&](const md::stack::Model::Ptr&) -> mfe::DpState { return md::stack::DpState{}; },
  };
  return std::visit(vis, m);
}

pfn::PfnState CreatePfnState(const BackendModelPtr& m) {
  auto vis = overloaded{
      [&](const md::base::Model::Ptr&) -> pfn::PfnState { return md::base::PfnState{}; },
      [&](const md::base::opt::Model::Ptr&) -> pfn::PfnState { return md::base::PfnState{}; },
      [&](const md::stack::Model::Ptr&) -> pfn::PfnState { fatal("unimplemented"); },
  };
  return std::visit(vis, m);
}

}  // namespace

Ctx::Ctx(Ctx&& o) noexcept
    : cfg_(std::move(o.cfg_)), backend_(o.backend_), backends_(std::move(o.backends_)),
      backend_once_{} {}

Ctx& Ctx::operator=(Ctx&& o) noexcept {
  if (this == &o) return *this;
  this->~Ctx();
  new (this) Ctx(std::move(o));
  return *this;
}

const BackendModelPtr& Ctx::EnsureBackend(BackendKind kind) const {
  auto idx = static_cast<size_t>(kind);
  std::call_once(backend_once_[idx], [this, idx, kind]() {
    if (backends_[idx].has_value()) return;
    backends_[idx] = BackendFromBackendCfg(kind, cfg_);
  });
  verify(backends_[idx].has_value(), "no backend available");
  return *backends_[idx];
}

MfeBackend Ctx::BackendForFold(
    std::optional<MfeAlg> alg, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf) const {
  if (alg.has_value() && *alg == MfeAlg::BRUTE) {
    std::string log;
    auto resolved = ResolveEfn(backend_, cfg_, cfg, pf, &log);
    if (!resolved.has_value())
      fatal("No backend supports this configuration for mfe brute:\n{}", log);
    const auto& m = EnsureBackend(resolved->backend);
    return {.m = m,
        .alg = MfeAlg::BRUTE,
        .mfe_fn = nullptr,
        .mfe_exterior_fn = nullptr,
        .trace_fn = nullptr};
  }

  std::string log;
  auto resolved = ResolveMfe(backend_, alg, cfg_, cfg, pf, &log);
  if (!resolved.has_value()) fatal("No MFE algorithm supports configuration:\n{}", log);

  const auto& m = EnsureBackend(resolved->backend);
  return {.m = m,
      .alg = resolved->alg,
      .mfe_fn = GetMfeFn(resolved->backend, resolved->alg),
      .mfe_exterior_fn = GetMfeExteriorFn(resolved->backend),
      .trace_fn = GetTraceFn(resolved->backend)};
}

SuboptBackend Ctx::BackendForSubopt(std::optional<SuboptAlg> alg, std::optional<MfeAlg> mfe_alg,
    const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf,
    const subopt::SuboptCfg& subopt_cfg) const {
  if (alg.has_value() && *alg == SuboptAlg::BRUTE) {
    std::string log;
    auto resolved = ResolveEfn(backend_, cfg_, cfg, pf, &log);
    if (!resolved.has_value())
      fatal("No backend supports this configuration for subopt brute:\n{}", log);
    const auto& m = EnsureBackend(resolved->backend);
    return {.m = m,
        .alg = SuboptAlg::BRUTE,
        .mfe_alg = mfe_alg.value_or(MfeAlg::BRUTE),
        .subopt_fn = nullptr,
        .mfe_fn = nullptr,
        .mfe_exterior_fn = nullptr};
  }

  auto mfe_backend = BackendForFold(mfe_alg, cfg, pf);
  auto kind = GetBackendKind(mfe_backend.m);

  std::string log;
  auto resolved = ResolveSubopt(kind, alg, cfg_, cfg, pf, subopt_cfg, &log);
  if (!resolved.has_value()) fatal("No subopt algorithm supports configuration:\n{}", log);

  if (mfe_backend.alg == MfeAlg::BRUTE) {
    fatal("subopt algorithm {} does not support mfe_alg={}", resolved->alg, mfe_backend.alg);
  }
  return {.m = mfe_backend.m,
      .alg = resolved->alg,
      .mfe_alg = mfe_backend.alg,
      .subopt_fn = GetSuboptFn(kind, resolved->alg),
      .mfe_fn = mfe_backend.mfe_fn,
      .mfe_exterior_fn = mfe_backend.mfe_exterior_fn};
}

PfnBackend Ctx::BackendForPfn(
    std::optional<PfnAlg> alg, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf) const {
  if (alg.has_value() && *alg == PfnAlg::BRUTE) {
    std::string log;
    auto resolved = ResolveEfn(backend_, cfg_, cfg, pf, &log);
    if (!resolved.has_value())
      fatal("No backend supports this configuration for pfn brute:\n{}", log);
    const auto& m = EnsureBackend(resolved->backend);
    return {.m = m, .alg = PfnAlg::BRUTE, .pfn_fn = nullptr};
  }

  std::string log;
  auto resolved = ResolvePfn(backend_, alg, cfg_, cfg, pf, &log);
  if (!resolved.has_value()) fatal("No PFN algorithm supports configuration:\n{}", log);

  const auto& m = EnsureBackend(resolved->backend);
  return {.m = m, .alg = resolved->alg, .pfn_fn = GetPfnFn(resolved->backend, resolved->alg)};
}

erg::EnergyResult Ctx::Efn(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, const Ctds* given_ctd, bool build_structure) const {
  std::string log;
  auto resolved = ResolveEfn(backend_, cfg_, cfg, pf, &log);
  if (!resolved.has_value()) fatal("No backend supports this configuration for Efn:\n{}", log);
  const auto& m = EnsureBackend(resolved->backend);
  return TotalEnergy(m, r, s, given_ctd, cfg, pf, build_structure);
}

FoldResult Ctx::Fold(const Primary& r, std::optional<MfeAlg> alg, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) const {
  auto backend = BackendForFold(alg, cfg, pf);

  if (backend.alg == MfeAlg::BRUTE) {
    auto subopt = md::brute::MfeBrute(r, backend.m, cfg_, cfg, pf);
    return {.mfe = {.dp{}, .energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  mfe::DpState dp = CreateDpState(backend.m);
  backend.mfe_fn(backend.m, r, dp, cfg, pf);
  auto energy = backend.mfe_exterior_fn(backend.m, r, dp, cfg, pf);
  auto tb = backend.trace_fn(backend.m, r, dp, cfg, pf, trace_cfg);
  return FoldResult{
      .mfe = {.dp = std::move(dp), .energy = energy},
      .tb = std::move(tb),
  };
}

std::vector<subopt::SuboptResult> Ctx::SuboptIntoVector(const Primary& r,
    std::optional<MfeAlg> mfe_alg, std::optional<SuboptAlg> alg, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, subopt::SuboptCfg subopt_cfg) const {
  std::vector<subopt::SuboptResult> subopts;
  [[maybe_unused]] const int strucs = Subopt(
      r, mfe_alg, alg, cfg, pf,
      [&subopts](const subopt::SuboptResult& subopt) { subopts.push_back(subopt); }, subopt_cfg);
  assert(strucs == static_cast<int>(subopts.size()));
  return subopts;
}

int Ctx::Subopt(const Primary& r, std::optional<MfeAlg> mfe_alg, std::optional<SuboptAlg> alg,
    erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
    subopt::SuboptCfg subopt_cfg) const {
  auto backend = BackendForSubopt(alg, mfe_alg, cfg, pf, subopt_cfg);

  if (backend.alg == SuboptAlg::BRUTE) {
    auto subopts = md::brute::SuboptBrute(r, backend.m, cfg_, cfg, pf, subopt_cfg);
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  mfe::DpState dp = CreateDpState(backend.m);
  backend.mfe_fn(backend.m, r, dp, cfg, pf);
  backend.mfe_exterior_fn(backend.m, r, dp, cfg, pf);

  return backend.subopt_fn(backend.m, Primary(r), std::move(dp), cfg, pf, fn, subopt_cfg);
}

pfn::PfnResult Ctx::Pfn(const Primary& r, std::optional<PfnAlg> alg, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf) const {
  auto backend = BackendForPfn(alg, cfg, pf);

  if (backend.alg == PfnAlg::BRUTE) return md::brute::PfnBrute(r, backend.m, cfg_, cfg, pf);

  pfn::PfnState dp = CreatePfnState(backend.m);
  auto pfn = backend.pfn_fn(backend.m, r, dp, cfg, pf);

  return pfn::PfnResult{.state = std::move(dp), .pfn = std::move(pfn)};
}

Ctx Ctx::FromArgParse(const ArgParse& args) {
  return Ctx(args.MaybeGet<BackendKind>(OPT_BACKEND), BackendCfg::FromArgParse(args));
}

void RegisterOpts(ArgParse* args) {
  RegisterOptsBackendCfg(args);
  RegisterOptsAlgorithm(args);
  erg::RegisterOptsPseudofree(args);
  trace::RegisterOpts(args);
  subopt::RegisterOpts(args);
}

}  // namespace mrna
