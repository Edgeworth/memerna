// Copyright 2016 Eliot Courtney.
#include "api/ctx/ctx.h"

#include <cassert>
#include <new>
#include <utility>
#include <variant>
#include <vector>

#include "api/ctx/backend.h"
#include "api/mfe.h"
#include "api/pfn.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace_cfg.h"
#include "backends/base/energy/boltz_model.h"
#include "backends/base/mfe/mfe_debug.h"
#include "backends/base/mfe/mfe_exterior.h"
#include "backends/base/mfe/mfe_lyngso_sparse_opt.h"
#include "backends/base/mfe/mfe_opt.h"
#include "backends/base/mfe/mfe_sparse_opt.h"
#include "backends/base/pfn/pfn_debug.h"
#include "backends/base/pfn/pfn_opt.h"
#include "backends/base/subopt/subopt_debug.h"
#include "backends/base/subopt/subopt_iterative.h"
#include "backends/base/subopt/subopt_persistent.h"
#include "backends/base/trace/trace.h"
#include "backends/baseopt/energy/boltz_model.h"
#include "backends/baseopt/mfe/mfe_debug.h"
#include "backends/baseopt/mfe/mfe_exterior.h"
#include "backends/baseopt/mfe/mfe_lyngso_sparse_opt.h"
#include "backends/baseopt/mfe/mfe_opt.h"
#include "backends/baseopt/mfe/mfe_sparse_opt.h"
#include "backends/baseopt/pfn/pfn_debug.h"
#include "backends/baseopt/pfn/pfn_opt.h"
#include "backends/baseopt/subopt/subopt_debug.h"
#include "backends/baseopt/subopt/subopt_iterative.h"
#include "backends/baseopt/subopt/subopt_persistent.h"
#include "backends/baseopt/trace/trace.h"
#include "backends/brute/alg.h"
#include "backends/common/base/dp.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/dp.h"
#include "backends/stack/mfe/mfe_opt.h"
#include "backends/stack/mfe/mfe_exterior.h"
#include "backends/stack/subopt/subopt_iterative.h"
#include "backends/stack/subopt/subopt_persistent.h"
#include "backends/stack/trace/trace.h"
#include "model/energy.h"
#include "model/pfn.h"
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
    : cfg_(std::move(o.cfg_)), backends_(std::move(o.backends_)), backend_once_{} {}

Ctx& Ctx::operator=(Ctx&& o) noexcept {
  if (this == &o) return *this;
  this->~Ctx();
  new (this) Ctx(std::move(o));
  return *this;
}

const BackendModelPtr& Ctx::EnsureBackend() const {
  if (!cfg_.has_value()) {
    // Injected model mode - model already at backends_[0]
    verify(backends_[0].has_value(), "no backend available");
    return *backends_[0];
  }
  // Normal lazy loading mode
  auto idx = static_cast<size_t>(cfg_->backend);
  std::call_once(backend_once_[idx], [this, idx]() {
    if (backends_[idx].has_value()) return;
    backends_[idx] = BackendFromBackendCfg(*cfg_);
  });
  verify(backends_[idx].has_value(), "no backend available");
  return *backends_[idx];
}

const BackendModelPtr& Ctx::BackendForFold([[maybe_unused]] MfeAlg alg) const {
  return EnsureBackend();
}

const BackendModelPtr& Ctx::BackendForSubopt(
    [[maybe_unused]] SuboptAlg alg, [[maybe_unused]] MfeAlg mfe_alg) const {
  return EnsureBackend();
}

const BackendModelPtr& Ctx::BackendForPfn([[maybe_unused]] PfnAlg alg) const {
  return EnsureBackend();
}

erg::EnergyResult Ctx::Efn(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, const Ctds* given_ctd, bool build_structure) const {
  const auto& m = EnsureBackend();
  return TotalEnergy(m, r, s, given_ctd, cfg, pf, build_structure);
}

void Ctx::ComputeMfe(const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, MfeAlg alg,
    erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf) const {
  auto vis = overloaded{
      [&](const md::base::Model::Ptr& m) {
        auto& state = std::get<md::base::DpState>(dp);
        switch (alg) {
        case MfeAlg::DEBUG: md::base::MfeDebug::Run(r, m, state, cfg, pf); break;
        case MfeAlg::OPT: md::base::MfeOpt::Run(r, m, state, cfg, pf); break;
        case MfeAlg::AUTO:
        case MfeAlg::SPARSE_OPT: md::base::MfeSparseOpt::Run(r, m, state, cfg, pf); break;
        case MfeAlg::LYNGSO_SPARSE_OPT:
          md::base::MfeLyngsoSparseOpt::Run(r, m, state, cfg, pf);
          break;
        default: fatal("unsupported mfe algorithm for energy model: {}", alg);
        }
      },
      [&](const md::base::opt::Model::Ptr& m) {
        auto& state = std::get<md::base::DpState>(dp);
        switch (alg) {
        case MfeAlg::DEBUG: md::base::opt::MfeDebug::Run(r, m, state, cfg, pf); break;
        case MfeAlg::OPT: md::base::opt::MfeOpt::Run(r, m, state, cfg, pf); break;
        case MfeAlg::AUTO:
        case MfeAlg::SPARSE_OPT: md::base::opt::MfeSparseOpt::Run(r, m, state, cfg, pf); break;
        case MfeAlg::LYNGSO_SPARSE_OPT:
          md::base::opt::MfeLyngsoSparseOpt::Run(r, m, state, cfg, pf);
          break;
        default: fatal("unsupported mfe algorithm for energy model: {}", alg);
        }
      },
      [&](const md::stack::Model::Ptr& m) {
        auto& state = std::get<md::stack::DpState>(dp);
        switch (alg) {
        case MfeAlg::AUTO:
        case MfeAlg::OPT: md::stack::MfeOpt::Run(r, m, state, cfg, pf); break;
        default: fatal("unsupported mfe algorithm for energy model: {}", alg);
        }
      },
  };
  std::visit(vis, m);
}

Energy Ctx::ComputeMfeExterior(const BackendModelPtr& m, const Primary& r, mfe::DpState& dp,
    erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf) const {
  auto vis = overloaded{
      [&](const md::base::Model::Ptr& m) {
        auto& state = std::get<md::base::DpState>(dp);
        return md::base::MfeExterior(r, m, state, cfg, pf);
      },
      [&](const md::base::opt::Model::Ptr& m) {
        auto& state = std::get<md::base::DpState>(dp);
        return md::base::opt::MfeExterior(r, m, state, cfg, pf);
      },
      [&](const md::stack::Model::Ptr& m) {
        auto& state = std::get<md::stack::DpState>(dp);
        return md::stack::MfeExterior(r, m, state, cfg, pf);
      },
  };
  return std::visit(vis, m);
}

trace::TraceResult Ctx::ComputeTraceback(const BackendModelPtr& m, const Primary& r,
    const mfe::DpState& dp, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
    const trace::TraceCfg& trace_cfg) const {
  auto vis = overloaded{
      [&](const md::base::Model::Ptr& m) -> trace::TraceResult {
        const auto& state = std::get<md::base::DpState>(dp);
        return md::base::Traceback(r, m, state, cfg, pf, trace_cfg);
      },
      [&](const md::base::opt::Model::Ptr& m) {
        const auto& state = std::get<md::base::DpState>(dp);
        return md::base::opt::Traceback(r, m, state, cfg, pf, trace_cfg);
      },
      [&](const md::stack::Model::Ptr& m) -> trace::TraceResult {
        const auto& state = std::get<md::stack::DpState>(dp);
        return md::stack::Traceback(r, m, state, cfg, pf, trace_cfg);
      },
  };
  return std::visit(vis, m);
}

FoldResult Ctx::Fold(const Primary& r, MfeAlg alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
    const trace::TraceCfg& trace_cfg) const {
  const auto& m = BackendForFold(alg);

  // Resolve AUTO to the default for this backend
  if (alg == MfeAlg::AUTO) {
    alg = (GetBackendKind(m) == BackendKind::STACK) ? MfeAlg::OPT : MfeAlg::SPARSE_OPT;
  }

  if (alg == MfeAlg::BRUTE) {
    auto subopt = md::brute::MfeBrute(r, m, cfg, pf);
    return {.mfe = {.dp{}, .energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  mfe::DpState dp = CreateDpState(m);
  ComputeMfe(m, r, dp, alg, cfg, pf);
  auto energy = ComputeMfeExterior(m, r, dp, cfg, pf);
  auto tb = ComputeTraceback(m, r, dp, cfg, pf, trace_cfg);
  return FoldResult{
      .mfe = {.dp = std::move(dp), .energy = energy},
      .tb = std::move(tb),
  };
}

std::vector<subopt::SuboptResult> Ctx::SuboptIntoVector(const Primary& r, MfeAlg mfe_alg,
    SuboptAlg alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
    subopt::SuboptCfg subopt_cfg) const {
  std::vector<subopt::SuboptResult> subopts;
  [[maybe_unused]] const int strucs = Subopt(
      r, mfe_alg, alg, cfg, pf,
      [&subopts](const subopt::SuboptResult& subopt) { subopts.push_back(subopt); }, subopt_cfg);
  assert(strucs == static_cast<int>(subopts.size()));
  return subopts;
}

int Ctx::Subopt(const Primary& r, MfeAlg mfe_alg, SuboptAlg alg, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
    subopt::SuboptCfg subopt_cfg) const {
  const auto& m = BackendForSubopt(alg, mfe_alg);

  // Resolve AUTO to the default for this backend
  if (alg == SuboptAlg::AUTO) {
    alg = SuboptAlg::ITERATIVE;
  }

  if (alg == SuboptAlg::BRUTE) {
    // TODO(3): handle cases other than max structures.
    auto subopts = md::brute::SuboptBrute(r, m, cfg, pf, subopt_cfg);
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  // Resolve MfeAlg::AUTO to the default for this backend
  if (mfe_alg == MfeAlg::AUTO) {
    mfe_alg = (GetBackendKind(m) == BackendKind::STACK) ? MfeAlg::OPT : MfeAlg::SPARSE_OPT;
  }

  mfe::DpState dp = CreateDpState(m);
  ComputeMfe(m, r, dp, mfe_alg, cfg, pf);
  ComputeMfeExterior(m, r, dp, cfg, pf);

  auto vis = overloaded{//
      [&](const md::base::Model::Ptr& m) mutable -> int {
        auto state = std::get<md::base::DpState>(std::move(dp));
        switch (alg) {
        case SuboptAlg::DEBUG:
          return md::base::SuboptDebug(Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::ITERATIVE:
          return md::base::SuboptIterative<false>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::ITERATIVE_LOWMEM:
          return md::base::SuboptIterative<true>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::PERSISTENT:
          return md::base::SuboptPersistent<false>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::PERSISTENT_LOWMEM:
          return md::base::SuboptPersistent<true>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        default: fatal("unsupported subopt algorithm for energy model: {}", alg);
        }
      },
      [&](const md::base::opt::Model::Ptr& m) mutable -> int {
        auto state = std::get<md::base::DpState>(std::move(dp));
        switch (alg) {
        case SuboptAlg::DEBUG:
          return md::base::opt::SuboptDebug(Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::ITERATIVE:
          return md::base::opt::SuboptIterative<false>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::ITERATIVE_LOWMEM:
          return md::base::opt::SuboptIterative<true>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::PERSISTENT:
          return md::base::opt::SuboptPersistent<false>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::PERSISTENT_LOWMEM:
          return md::base::opt::SuboptPersistent<true>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        default: fatal("unsupported subopt algorithm for energy model: {}", alg);
        }
      },
      [&](const md::stack::Model::Ptr& m) mutable -> int {
        auto state = std::get<md::stack::DpState>(std::move(dp));
        switch (alg) {
        case SuboptAlg::ITERATIVE:
          return md::stack::SuboptIterative<false>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::ITERATIVE_LOWMEM:
          return md::stack::SuboptIterative<true>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::PERSISTENT:
          return md::stack::SuboptPersistent<false>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        case SuboptAlg::PERSISTENT_LOWMEM:
          return md::stack::SuboptPersistent<true>(
              Primary(r), m, std::move(state), cfg, pf, subopt_cfg)
              .Run(fn);
        default: fatal("unsupported subopt algorithm for energy model: {}", alg);
        }
      }};
  return std::visit(vis, m);
}

pfn::PfnResult Ctx::Pfn(
    const Primary& r, PfnAlg alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf) const {
  const auto& m = BackendForPfn(alg);

  // PFN algorithms require bulge_states = false
  cfg.bulge_states = false;

  // Resolve AUTO to the default for this backend
  if (alg == PfnAlg::AUTO) {
    alg = (GetBackendKind(m) == BackendKind::STACK) ? PfnAlg::BRUTE : PfnAlg::OPT;
  }

  if (alg == PfnAlg::BRUTE) {
    return md::brute::PfnBrute(r, m, cfg, pf);
  }

  pfn::PfnState dp = CreatePfnState(m);

  auto vis = overloaded{
      [&](const md::base::Model::Ptr& m) -> PfnTables {
        auto state = std::get<md::base::PfnState>(std::move(dp));
        switch (alg) {
        case PfnAlg::DEBUG: return md::base::PfnDebug::Run(r, m, cfg, state, pf);
        case PfnAlg::OPT:
          return md::base::PfnOpt::Run(r, md::base::BoltzModel::Create(m), cfg, state, pf);
        default: fatal("unsupported partition algorithm for energy model: {}", alg);
        }
      },
      [&](const md::base::opt::Model::Ptr& m) -> PfnTables {
        auto state = std::get<md::base::PfnState>(std::move(dp));
        switch (alg) {
        case PfnAlg::DEBUG: return md::base::opt::PfnDebug::Run(r, m, cfg, state, pf);
        case PfnAlg::OPT:
          return md::base::opt::PfnOpt::Run(
              r, md::base::opt::BoltzModel::Create(m), cfg, state, pf);
        default: fatal("unsupported partition algorithm for energy model: {}", alg);
        }
      },
      // TODO(2): Implement partition for t22.
      [&](const md::stack::Model::Ptr&) -> PfnTables { fatal("unimplemented"); },
  };
  auto pfn = std::visit(vis, m);

  return pfn::PfnResult{.state = std::move(dp), .pfn = std::move(pfn)};
}

Ctx Ctx::FromArgParse(const ArgParse& args) { return Ctx(BackendCfg::FromArgParse(args)); }

void RegisterOpts(ArgParse* args) {
  RegisterOptsBackendCfg(args);
  RegisterOptsAlgorithm(args);
  erg::RegisterOptsPseudofree(args);
  trace::RegisterOpts(args);
  subopt::RegisterOpts(args);
}

}  // namespace mrna
