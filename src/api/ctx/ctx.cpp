// Copyright 2016 Eliot Courtney.
#include "api/ctx/ctx.h"

#include <algorithm>
#include <cassert>
#include <utility>
#include <variant>
#include <vector>

#include "api/ctx/ctx_cfg.h"
#include "api/mfe.h"
#include "api/part.h"
#include "api/trace/trace_cfg.h"
#include "model/energy.h"
#include "model/part.h"
#include "model/primary.h"
#include "models/brute/alg.h"
#include "models/t04/energy/boltz_model.h"
#include "models/t04/energy/model.h"
#include "models/t04/mfe/dp.h"
#include "models/t04/mfe/mfe.h"
#include "models/t04/part/part.h"
#include "models/t04/subopt/subopt_debug.h"
#include "models/t04/subopt/subopt_iterative.h"
#include "models/t04/subopt/subopt_persistent.h"
#include "models/t04/trace/trace.h"
#include "models/t22/energy/model.h"
#include "models/t22/mfe/mfe.h"
#include "models/t22/subopt/subopt_iterative.h"
#include "models/t22/subopt/subopt_persistent.h"
#include "models/t22/trace/trace.h"
#include "util/error.h"
#include "util/util.h"

namespace mrna {

namespace {

mfe::DpState CreateDpState(const erg::EnergyModelPtr& em) {
  auto vis = overloaded{
      [&](const md::t04::Model::Ptr&) -> mfe::DpState { return md::t04::DpState{}; },
      [&](const md::t22::Model::Ptr&) -> mfe::DpState { return md::t22::DpState{}; },
  };
  return std::visit(vis, em);
}

part::PartState CreatePartState(const erg::EnergyModelPtr& em) {
  auto vis = overloaded{
      [&](const md::t04::Model::Ptr&) -> part::PartState { return md::t04::PartState{}; },
      [&](const md::t22::Model::Ptr&) -> part::PartState { fatal("unimplemented"); },
  };
  return std::visit(vis, em);
}

}  // namespace

erg::EnergyResult Ctx::Efn(
    const Primary& r, const Secondary& s, const Ctds* given_ctd, bool build_structure) const {
  return erg::TotalEnergy(em(), r, s, given_ctd, build_structure);
}

void Ctx::ComputeMfe(const Primary& r, mfe::DpState& dp) const {
  auto vis = overloaded{
      [&](const md::t04::Model::Ptr& em) {
        auto& state = std::get<md::t04::DpState>(dp);
        switch (cfg_.mfe_alg) {
        case CtxCfg::MfeAlg::DEBUG: md::t04::MfeDebug(r, em, state); break;
        case CtxCfg::MfeAlg::OPT: md::t04::MfeOpt(r, em, state); break;
        case CtxCfg::MfeAlg::AUTO:
        case CtxCfg::MfeAlg::SPARSE_OPT: md::t04::MfeSparseOpt(r, em, state); break;
        case CtxCfg::MfeAlg::LYNGSO_SPARSE_OPT: md::t04::MfeLyngsoSparseOpt(r, em, state); break;
        default: fatal("unsupported mfe algorithm for energy model: {}", cfg_.mfe_alg);
        }
      },
      [&](const md::t22::Model::Ptr& em) {
        auto& state = std::get<md::t22::DpState>(dp);
        switch (cfg_.mfe_alg) {
        case CtxCfg::MfeAlg::AUTO:
        case CtxCfg::MfeAlg::DEBUG: md::t22::MfeDebug(r, em, state); break;
        default: fatal("unsupported mfe algorithm for energy model: {}", cfg_.mfe_alg);
        }
      },
  };
  std::visit(vis, em_);
}

Energy Ctx::ComputeMfeExterior(const Primary& r, mfe::DpState& dp) const {
  auto vis = overloaded{
      [&](const md::t04::Model::Ptr& em) {
        auto& state = std::get<md::t04::DpState>(dp);
        return md::t04::MfeExterior(r, em, state);
      },
      [&](const md::t22::Model::Ptr& em) {
        auto& state = std::get<md::t22::DpState>(dp);
        return md::t22::MfeExterior(r, em, state);
      },
  };
  return std::visit(vis, em_);
}

trace::TraceResult Ctx::ComputeTraceback(
    const Primary& r, const trace::TraceCfg& cfg, const mfe::DpState& dp) const {
  auto vis = overloaded{
      [&](const md::t04::Model::Ptr& em) -> trace::TraceResult {
        const auto& state = std::get<md::t04::DpState>(dp);
        return md::t04::Traceback(r, em, cfg, state);
      },
      [&](const md::t22::Model::Ptr& em) -> trace::TraceResult {
        const auto& state = std::get<md::t22::DpState>(dp);
        return md::t22::Traceback(r, em, cfg, state);
      },
  };
  return std::visit(vis, em_);
}

FoldResult Ctx::Fold(const Primary& r, const trace::TraceCfg& cfg) const {
  if (cfg_.mfe_alg == CtxCfg::MfeAlg::BRUTE) {
    auto subopt = md::brute::MfeBrute(r, em_);
    return {.mfe = {.dp{}, .energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  mfe::DpState dp = CreateDpState(em_);
  ComputeMfe(r, dp);
  auto energy = ComputeMfeExterior(r, dp);
  auto tb = ComputeTraceback(r, cfg, dp);
  return FoldResult{
      .mfe = {.dp = std::move(dp), .energy = energy},
      .tb = std::move(tb),
  };
}

std::vector<subopt::SuboptResult> Ctx::SuboptimalIntoVector(
    const Primary& r, subopt::SuboptCfg cfg) const {
  std::vector<subopt::SuboptResult> subopts;
  [[maybe_unused]] const int strucs = Suboptimal(
      r, [&subopts](const subopt::SuboptResult& subopt) { subopts.push_back(subopt); }, cfg);
  assert(strucs == static_cast<int>(subopts.size()));
  return subopts;
}

int Ctx::Suboptimal(
    const Primary& r, const subopt::SuboptCallback& fn, subopt::SuboptCfg cfg) const {
  if (cfg_.subopt_alg == CtxCfg::SuboptAlg::BRUTE) {
    // TODO(3): handle cases other than max structures.
    auto subopts = md::brute::SuboptBrute(r, em_, cfg);
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  mfe::DpState dp = CreateDpState(em_);
  ComputeMfe(r, dp);
  ComputeMfeExterior(r, dp);

  auto vis = overloaded{// fix clang-format comment
      [&](const md::t04::Model::Ptr& em) mutable -> int {
        auto state = std::get<md::t04::DpState>(std::move(dp));
        switch (cfg_.subopt_alg) {
        case CtxCfg::SuboptAlg::DEBUG:
          return md::t04::SuboptDebug(Primary(r), em, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::AUTO:
        case CtxCfg::SuboptAlg::ITERATIVE:
          return md::t04::SuboptIterative(Primary(r), em, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::PERSISTENT:
          return md::t04::SuboptPersistent(Primary(r), em, std::move(state), cfg).Run(fn);
        default: fatal("unsupported subopt algorithm for energy model: {}", cfg_.subopt_alg);
        }
      },
      [&](const md::t22::Model::Ptr& em) mutable -> int {
        auto state = std::get<md::t22::DpState>(std::move(dp));
        switch (cfg_.subopt_alg) {
        case CtxCfg::SuboptAlg::AUTO:
        case CtxCfg::SuboptAlg::ITERATIVE:
          return md::t22::SuboptIterative(Primary(r), em, std::move(state), cfg).Run(fn);
        case CtxCfg::SuboptAlg::PERSISTENT:
          return md::t22::SuboptPersistent(Primary(r), em, std::move(state), cfg).Run(fn);
        default: fatal("unsupported subopt algorithm for energy model: {}", cfg_.subopt_alg);
        }
      }};
  return std::visit(vis, em_);
}

part::PartResult Ctx::Partition(const Primary& r) const {
  if (cfg_.part_alg == CtxCfg::PartAlg::BRUTE) {
    return md::brute::PartitionBrute(r, em_);
  }

  part::PartState dp = CreatePartState(em_);

  auto vis = overloaded{
      [&](const md::t04::Model::Ptr& em) -> Part {
        auto state = std::get<md::t04::PartState>(std::move(dp));
        switch (cfg_.part_alg) {
        case CtxCfg::PartAlg::DEBUG: return md::t04::PartitionDebug(r, em, state);
        case CtxCfg::PartAlg::AUTO:
        case CtxCfg::PartAlg::OPT:
          return md::t04::PartitionOpt(r, md::t04::BoltzModel::Create(em), state);
        default: fatal("unsupported partition algorithm for energy model: {}", cfg_.part_alg);
        }
      },
      // TODO(2): Implement partition for t22.
      [&](const md::t22::Model::Ptr&) -> Part { fatal("unimplemented"); },
  };
  auto part = std::visit(vis, em_);

  return part::PartResult{.state = std::move(dp), .part = std::move(part)};
}

Ctx Ctx::FromArgParse(const ArgParse& args) {
  return {erg::FromArgParse(args), CtxCfg::FromArgParse(args)};
}

}  // namespace mrna
