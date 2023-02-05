// Copyright 2016 Eliot Courtney.
#include "api/ctx/ctx.h"

#include <cassert>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "api/part.h"
#include "model/energy.h"
#include "model/primary.h"
#include "models/brute/alg.h"
#include "models/t04/energy/boltz_model.h"
#include "models/t04/energy/model.h"
#include "models/t04/mfe/dp.h"
#include "models/t04/mfe/mfe.h"
#include "models/t04/part/part.h"
#include "models/t04/part/partition.h"
#include "models/t04/subopt/subopt_fastest.h"
#include "models/t04/subopt/subopt_slowest.h"
#include "models/t04/trace/trace.h"
#include "models/t22/energy/model.h"
#include "models/t22/mfe/mfe.h"
#include "models/t22/trace/trace.h"
#include "util/array.h"
#include "util/error.h"
#include "util/util.h"

namespace mrna::api {

erg::EnergyResult Ctx::Efn(
    const Primary& r, const Secondary& s, const Ctds* given_ctd, bool build_structure) const {
  return erg::TotalEnergy(em(), r, s, given_ctd, build_structure);
}

DpArray Ctx::ComputeMfe(const Primary& r) const {
  auto vis = overloaded{
      [&](const erg::t04::Model::Ptr& em) -> DpArray {
        switch (cfg_.dp_alg) {
        case CtxCfg::DpAlg::SLOWEST: return mfe::t04::MfeSlowest(r, em);
        case CtxCfg::DpAlg::SLOW: return mfe::t04::MfeSlow(r, em);
        case CtxCfg::DpAlg::FASTEST: return mfe::t04::MfeFastest(r, em);
        case CtxCfg::DpAlg::LYNGSO: return mfe::t04::MfeLyngso(r, em);
        default: bug();
        }
      },
      [&](const erg::t22::Model::Ptr& em) {
        switch (cfg_.dp_alg) {
        case CtxCfg::DpAlg::SLOWEST:
        case CtxCfg::DpAlg::SLOW:
        case CtxCfg::DpAlg::FASTEST:
        case CtxCfg::DpAlg::LYNGSO: return mfe::t22::MfeSlowest(r, em);
        default: bug();
        }
      },
  };

  return std::visit(vis,

      em_);
}

ExtArray Ctx::ComputeMfeExterior(const Primary& r, const DpArray& dp) const {
  auto vis = overloaded{
      [&](const erg::t04::Model::Ptr& em) -> ExtArray { return mfe::t04::MfeExterior(r, em, dp); },
      [&](const erg::t22::Model::Ptr& em) -> ExtArray { return mfe::t22::MfeExterior(r, em, dp); },
  };
  return std::visit(vis, em_);
}

trace::TraceResult Ctx::ComputeTraceback(
    const Primary& r, const DpArray& dp, const ExtArray& ext) const {
  auto vis = overloaded{
      [&](const erg::t04::Model::Ptr& em) -> trace::TraceResult {
        return tb::t04::Traceback(r, em, dp, ext);
      },
      [&](const erg::t22::Model::Ptr& em) -> trace::TraceResult {
        return tb::t22::Traceback(r, em, dp, ext);
      },
  };
  return std::visit(vis, em_);
}

ctx::FoldResult Ctx::Fold(const Primary& r) const {
  if (cfg_.dp_alg == CtxCfg::DpAlg::BRUTE) {
    auto subopt = brute::MfeBrute(r, em_);
    return {.mfe = {.dp{}, .ext{}, .energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  auto dp = ComputeMfe(r);
  auto ext = ComputeMfeExterior(r, dp);
  auto energy = ext[0][EXT];
  auto tb = ComputeTraceback(r, dp, ext);
  return ctx::FoldResult{
      .mfe = {.dp = std::move(dp), .ext = std::move(ext), .energy = energy},
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
    auto subopts = brute::SuboptBrute(r, em_, cfg);
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  auto dp = ComputeMfe(r);
  auto ext = ComputeMfeExterior(r, dp);

  auto vis = overloaded{
      [&, ext = std::move(ext)](const erg::t04::Model::Ptr& em) mutable -> int {
        switch (cfg_.subopt_alg) {
        case CtxCfg::SuboptAlg::SLOWEST:
          return subopt::t04::SuboptSlowest(Primary(r), em, std::move(dp), std::move(ext), cfg)
              .Run(fn);
        case CtxCfg::SuboptAlg::FASTEST:
          return subopt::t04::SuboptFastest(Primary(r), em, std::move(dp), std::move(ext), cfg)
              .Run(fn);
        default: bug();
        }
      },
      // TODO(2): Implement suboptimal for t22.
      [&, ext = std::move(ext)](const erg::t22::Model::Ptr&) mutable -> int { bug(); },
  };
  return std::visit(vis, em_);
}

part::PartResult Ctx::Partition(const Primary& r) const {
  if (cfg_.part_alg == CtxCfg::PartAlg::BRUTE) {
    return brute::PartitionBrute(r, em_);
  }

  std::tuple<BoltzDpArray, BoltzExtArray> res;

  auto vis = overloaded{
      [&](const erg::t04::Model::Ptr& em) {
        switch (cfg_.part_alg) {
        case CtxCfg::PartAlg::SLOWEST: res = part::t04::PartitionSlowest(r, em); break;
        case CtxCfg::PartAlg::FASTEST:
          res = part::t04::PartitionFastest(r, erg::t04::BoltzModel::Create(em));
          break;
        default: bug();
        }
      },
      // TODO(2): Implement partition for t22.
      [&](const erg::t22::Model::Ptr&) { bug(); },
  };
  std::visit(vis, em_);

  const int N = static_cast<int>(r.size());
  auto [dp, ext] = std::move(res);
  BoltzSums p(N, 0);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) p[i][j] = dp[i][j][PT_P];

  part::Part part{std::move(p), ext[0][PTEXT_R]};
  auto prob = part.Prob();
  return part::PartResult{
      .dp = std::move(dp), .ext = std::move(ext), .part = std::move(part), .prob = std::move(prob)};
}

Ctx Ctx::FromArgParse(const ArgParse& args) {
  return {erg::FromArgParse(args), CtxCfg::FromArgParse(args)};
}

}  // namespace mrna::api
