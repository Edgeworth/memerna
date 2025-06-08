// Copyright 2016 Eliot Courtney.
#include "api/ctx/ctx.h"

#include <cassert>
#include <utility>
#include <variant>
#include <vector>

#include "api/ctx/backend.h"
#include "api/ctx/ctx_cfg.h"
#include "api/mfe.h"
#include "api/pfn.h"
#include "api/trace/trace_cfg.h"
#include "backends/base/energy/boltz_model.h"
#include "backends/base/mfe/mfe.h"
#include "backends/base/pfn/pfn.h"
#include "backends/base/subopt/subopt_debug.h"
#include "backends/base/subopt/subopt_iterative.h"
#include "backends/base/subopt/subopt_persistent.h"
#include "backends/base/trace/trace.h"
#include "backends/baseopt/energy/boltz_model.h"
#include "backends/baseopt/mfe/mfe.h"
#include "backends/baseopt/pfn/pfn.h"
#include "backends/baseopt/subopt/subopt_debug.h"
#include "backends/baseopt/subopt/subopt_iterative.h"
#include "backends/baseopt/subopt/subopt_persistent.h"
#include "backends/baseopt/trace/trace.h"
#include "backends/brute/alg.h"
#include "backends/common/base/dp.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/mfe.h"
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

erg::EnergyResult Ctx::Efn(
    const Primary& r, const Secondary& s, const Ctds* given_ctd, bool build_structure) const {
  return TotalEnergy(m(), r, s, given_ctd, build_structure);
}

void Ctx::ComputeMfe(const Primary& r, mfe::DpState& dp) const {
  auto vis = overloaded{
      [&](const md::base::Model::Ptr& m) {
        auto& state = std::get<md::base::DpState>(dp);
        switch (cfg_.mfe_alg) {
        case CtxCfg::MfeAlg::DEBUG: md::base::MfeDebug(r, m, state); break;
        case CtxCfg::MfeAlg::OPT: md::base::MfeOpt(r, m, state); break;
        case CtxCfg::MfeAlg::AUTO:
        case CtxCfg::MfeAlg::SPARSE_OPT: md::base::MfeSparseOpt(r, m, state); break;
        case CtxCfg::MfeAlg::LYNGSO_SPARSE_OPT: md::base::MfeLyngsoSparseOpt(r, m, state); break;
        default: fatal("unsupported mfe algorithm for energy model: {}", cfg_.mfe_alg);
        }
      },
      [&](const md::base::opt::Model::Ptr& m) {
        auto& state = std::get<md::base::DpState>(dp);
        switch (cfg_.mfe_alg) {
        case CtxCfg::MfeAlg::DEBUG: md::base::opt::MfeDebug(r, m, state); break;
        case CtxCfg::MfeAlg::OPT: md::base::opt::MfeOpt(r, m, state); break;
        case CtxCfg::MfeAlg::AUTO:
        case CtxCfg::MfeAlg::SPARSE_OPT: md::base::opt::MfeSparseOpt(r, m, state); break;
        case CtxCfg::MfeAlg::LYNGSO_SPARSE_OPT:
          md::base::opt::MfeLyngsoSparseOpt(r, m, state);
          break;
        default: fatal("unsupported mfe algorithm for energy model: {}", cfg_.mfe_alg);
        }
      },
      [&](const md::stack::Model::Ptr& m) {
        auto& state = std::get<md::stack::DpState>(dp);
        switch (cfg_.mfe_alg) {
        case CtxCfg::MfeAlg::AUTO:
        case CtxCfg::MfeAlg::DEBUG: md::stack::MfeDebug(r, m, state); break;
        default: fatal("unsupported mfe algorithm for energy model: {}", cfg_.mfe_alg);
        }
      },
  };
  std::visit(vis, m_);
}

Energy Ctx::ComputeMfeExterior(const Primary& r, mfe::DpState& dp) const {
  auto vis = overloaded{
      [&](const md::base::Model::Ptr& m) {
        auto& state = std::get<md::base::DpState>(dp);
        return md::base::MfeExterior(r, m, state);
      },
      [&](const md::base::opt::Model::Ptr& m) {
        auto& state = std::get<md::base::DpState>(dp);
        return md::base::opt::MfeExterior(r, m, state);
      },
      [&](const md::stack::Model::Ptr& m) {
        auto& state = std::get<md::stack::DpState>(dp);
        return md::stack::MfeExterior(r, m, state);
      },
  };
  return std::visit(vis, m_);
}

trace::TraceResult Ctx::ComputeTraceback(
    const Primary& r, const trace::TraceCfg& cfg, const mfe::DpState& dp) const {
  auto vis = overloaded{
      [&](const md::base::Model::Ptr& m) -> trace::TraceResult {
        const auto& state = std::get<md::base::DpState>(dp);
        return md::base::Traceback(r, m, cfg, state);
      },
      [&](const md::base::opt::Model::Ptr& m) {
        const auto& state = std::get<md::base::DpState>(dp);
        return md::base::opt::Traceback(r, m, cfg, state);
      },
      [&](const md::stack::Model::Ptr& m) -> trace::TraceResult {
        const auto& state = std::get<md::stack::DpState>(dp);
        return md::stack::Traceback(r, m, cfg, state);
      },
  };
  return std::visit(vis, m_);
}

FoldResult Ctx::Fold(const Primary& r, const trace::TraceCfg& cfg) const {
  if (cfg_.mfe_alg == CtxCfg::MfeAlg::BRUTE) {
    auto subopt = md::brute::MfeBrute(r, m_);
    return {.mfe = {.dp{}, .energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  mfe::DpState dp = CreateDpState(m_);
  ComputeMfe(r, dp);
  auto energy = ComputeMfeExterior(r, dp);
  auto tb = ComputeTraceback(r, cfg, dp);
  return FoldResult{
      .mfe = {.dp = std::move(dp), .energy = energy},
      .tb = std::move(tb),
  };
}

std::vector<subopt::SuboptResult> Ctx::SuboptIntoVector(
    const Primary& r, subopt::SuboptCfg cfg) const {
  std::vector<subopt::SuboptResult> subopts;
  [[maybe_unused]] const int strucs =
      Subopt(r, [&subopts](const subopt::SuboptResult& subopt) { subopts.push_back(subopt); }, cfg);
  assert(strucs == static_cast<int>(subopts.size()));
  return subopts;
}

int Ctx::Subopt(const Primary& r, const subopt::SuboptCallback& fn, subopt::SuboptCfg cfg) const {
  if (cfg_.subopt_alg == CtxCfg::SuboptAlg::BRUTE) {
    // TODO(3): handle cases other than max structures.
    auto subopts = md::brute::SuboptBrute(r, m_, cfg);
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  mfe::DpState dp = CreateDpState(m_);
  ComputeMfe(r, dp);
  ComputeMfeExterior(r, dp);

  auto vis = overloaded{
      [&](const md::base::Model::Ptr& m) mutable -> int {
        auto state = std::get<md::base::DpState>(std::move(dp));
        switch (cfg_.subopt_alg) {
        case CtxCfg::SuboptAlg::DEBUG:
          return md::base::SuboptDebug(Primary(r), m, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::AUTO:
        case CtxCfg::SuboptAlg::ITERATIVE:
          return md::base::SuboptIterative<false>(Primary(r), m, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::ITERATIVE_LOWMEM:
          return md::base::SuboptIterative<true>(Primary(r), m, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::PERSISTENT:
          return md::base::SuboptPersistent<false>(Primary(r), m, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::PERSISTENT_LOWMEM:
          return md::base::SuboptPersistent<true>(Primary(r), m, std::move(state), cfg).Run(fn);
        default: fatal("unsupported subopt algorithm for energy model: {}", cfg_.subopt_alg);
        }
      },
      [&](const md::base::opt::Model::Ptr& m) mutable -> int {
        auto state = std::get<md::base::DpState>(std::move(dp));
        switch (cfg_.subopt_alg) {
        case CtxCfg::SuboptAlg::DEBUG:
          return md::base::opt::SuboptDebug(Primary(r), m, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::AUTO:
        case CtxCfg::SuboptAlg::ITERATIVE:
          return md::base::opt::SuboptIterative<false>(Primary(r), m, std::move(state), cfg)
              .Run(fn);
        case CtxCfg::SuboptAlg::ITERATIVE_LOWMEM:
          return md::base::opt::SuboptIterative<true>(Primary(r), m, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::PERSISTENT:
          return md::base::opt::SuboptPersistent<true>(Primary(r), m, std::move(state), cfg)
              .Run(fn);
        case CtxCfg::SuboptAlg::PERSISTENT_LOWMEM:
          return md::base::opt::SuboptPersistent<true>(Primary(r), m, std::move(state), cfg)
              .Run(fn);
        default: fatal("unsupported subopt algorithm for energy model: {}", cfg_.subopt_alg);
        }
      },
      [&](const md::stack::Model::Ptr& m) mutable -> int {
        auto state = std::get<md::stack::DpState>(std::move(dp));
        switch (cfg_.subopt_alg) {
        case CtxCfg::SuboptAlg::AUTO:
        case CtxCfg::SuboptAlg::ITERATIVE:
          return md::stack::SuboptIterative<false>(Primary(r), m, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::ITERATIVE_LOWMEM:
          return md::stack::SuboptIterative<true>(Primary(r), m, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::PERSISTENT:
          return md::stack::SuboptPersistent<true>(Primary(r), m, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::PERSISTENT_LOWMEM:
          return md::stack::SuboptPersistent<true>(Primary(r), m, std::move(state), cfg).Run(fn);
        default: fatal("unsupported subopt algorithm for energy model: {}", cfg_.subopt_alg);
        }
      }};
  return std::visit(vis, m_);
}

pfn::PfnResult Ctx::Pfn(const Primary& r) const {
  if (cfg_.pfn_alg == CtxCfg::PfnAlg::BRUTE) {
    return md::brute::PfnBrute(r, m_);
  }

  pfn::PfnState dp = CreatePfnState(m_);

  auto vis = overloaded{
      [&](const md::base::Model::Ptr& m) -> PfnTables {
        auto state = std::get<md::base::PfnState>(std::move(dp));
        switch (cfg_.pfn_alg) {
        case CtxCfg::PfnAlg::DEBUG: return md::base::PfnDebug(r, m, state);
        case CtxCfg::PfnAlg::AUTO:
        case CtxCfg::PfnAlg::OPT:
          return md::base::PfnOpt(r, md::base::BoltzModel::Create(m), state);
        default: fatal("unsupported partition algorithm for energy model: {}", cfg_.pfn_alg);
        }
      },
      [&](const md::base::opt::Model::Ptr& m) -> PfnTables {
        auto state = std::get<md::base::PfnState>(std::move(dp));
        switch (cfg_.pfn_alg) {
        case CtxCfg::PfnAlg::DEBUG: return md::base::opt::PfnDebug(r, m, state);
        case CtxCfg::PfnAlg::AUTO:
        case CtxCfg::PfnAlg::OPT:
          return md::base::opt::PfnOpt(r, md::base::opt::BoltzModel::Create(m), state);
        default: fatal("unsupported partition algorithm for energy model: {}", cfg_.pfn_alg);
        }
      },
      // TODO(2): Implement partition for t22.
      [&](const md::stack::Model::Ptr&) -> PfnTables { fatal("unimplemented"); },
  };
  auto pfn = std::visit(vis, m_);

  return pfn::PfnResult{.state = std::move(dp), .pfn = std::move(pfn)};
}

Ctx Ctx::FromArgParse(const ArgParse& args) {
  return {BackendFromArgParse(args), CtxCfg::FromArgParse(args)};
}

}  // namespace mrna
