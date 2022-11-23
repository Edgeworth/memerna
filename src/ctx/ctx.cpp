// Copyright 2016 Eliot Courtney.
#include "ctx/ctx.h"

#include <cassert>
#include <tuple>
#include <utility>
#include <vector>

#include "compute/boltz_dp.h"
#include "compute/brute/alg.h"
#include "compute/dp.h"
#include "compute/energy/t04/boltz_model.h"
#include "compute/energy/t04/model.h"
#include "compute/mfe/mfe.h"
#include "compute/mfe/t04/mfe.h"
#include "compute/partition/partition.h"
#include "compute/subopt/t04/subopt_fastest.h"
#include "compute/subopt/t04/subopt_slowest.h"
#include "compute/traceback/t04/traceback.h"
#include "model/primary.h"
#include "util/array.h"
#include "util/error.h"

namespace mrna::ctx {

energy::EnergyResult Ctx::Efn(
    const Primary& r, const Secondary& s, const Ctds* given_ctd, bool build_structure) const {
  return energy::TotalEnergy(em(), r, s, given_ctd, build_structure);
}

DpArray Ctx::ComputeTables(const Primary& r) const {
  switch (energy::Kind(em())) {
  case energy::ModelKind::T04_LIKE:
    switch (cfg_.dp_alg) {
    case CtxCfg::DpAlg::SLOWEST: return mfe::t04::ComputeTablesSlowest(r, em_);
    case CtxCfg::DpAlg::SLOW: return mfe::t04::ComputeTablesSlow(r, em_);
    case CtxCfg::DpAlg::FASTEST: return mfe::t04::ComputeTablesFastest(r, em_);
    case CtxCfg::DpAlg::LYNGSO: return mfe::t04::ComputeTablesLyngso(r, em_);
    default: bug();
    }
    break;
  default: bug();
  }
}

std::tuple<ExtArray, Energy> Ctx::ComputeExterior(const Primary& r, const DpArray& dp) const {
  switch (energy::Kind(em())) {
  case energy::ModelKind::T04_LIKE:
    auto ext = mfe::t04::ComputeExterior(r, em(), dp);
    return {ext, ext[0][EXT]};
  default: bug();
  }
}

tb::TracebackResult Ctx::ComputeTraceback(
    const Primary& r, const DpArray& dp, const ExtArray& ext) const {
  switch (energy::Kind(em())) {
  case energy::ModelKind::T04_LIKE: return tb::t04::Traceback(r, em(), dp, ext); break;
  default: bug();
  }
}

ctx::FoldResult Ctx::Fold(const Primary& r) const {
  if (cfg_.dp_alg == CtxCfg::DpAlg::BRUTE) {
    auto subopt = brute::MfeBruteForce(r, em_);
    return {.mfe = {.dp{}, .ext{}, .energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  auto dp = ComputeTables(r);
  auto [ext, energy] = ComputeExterior(r, dp);
  auto tb = ComputeTraceback(r, dp, ext);
  return ctx::FoldResult{
      .mfe = {.dp = std::move(dp), .ext = std::move(ext), .energy = energy},
      .tb = std::move(tb),
  };
}

std::vector<subopt::SuboptResult> Ctx::SuboptimalIntoVector(
    const Primary& r, subopt::SuboptCfg cfg) const {
  std::vector<subopt::SuboptResult> subopts;
  [[maybe_unused]] int strucs = Suboptimal(
      r, [&subopts](const subopt::SuboptResult& subopt) { subopts.push_back(subopt); }, cfg);
  assert(strucs == static_cast<int>(subopts.size()));
  return subopts;
}

int Ctx::Suboptimal(
    const Primary& r, const subopt::SuboptCallback& fn, subopt::SuboptCfg cfg) const {
  if (cfg_.subopt_alg == CtxCfg::SuboptAlg::BRUTE) {
    // TODO(3): handle cases other than max structures.
    auto subopts = brute::SuboptimalBruteForce(r, em_, cfg);
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  auto dp = ComputeTables(r);
  auto [ext, energy] = ComputeExterior(r, dp);

  switch (energy::Kind(em())) {
  case energy::ModelKind::T04_LIKE:
    switch (cfg_.subopt_alg) {
    case CtxCfg::SuboptAlg::SLOWEST:
      return subopt::t04::SuboptimalSlowest(Primary(r), em_, std::move(dp), std::move(ext), cfg)
          .Run(fn);
    case CtxCfg::SuboptAlg::FASTEST:
      return subopt::t04::SuboptimalFastest(Primary(r), em_, std::move(dp), std::move(ext), cfg)
          .Run(fn);
    default: bug();
    }
    break;
  default: bug();
  }
}

part::PartResult Ctx::Partition(const Primary& r) const {
  if (cfg_.dp_alg == CtxCfg::DpAlg::BRUTE) {
    return part::PartitionBruteForce(r, em_);
  }

  std::tuple<BoltzDpArray, BoltzExtArray> res;
  switch (energy::Kind(em())) {
  case energy::ModelKind::T04_LIKE:
    switch (cfg_.part_alg) {
    case CtxCfg::PartAlg::SLOWEST: res = part::PartitionSlowest(r, em_); break;
    case CtxCfg::PartAlg::FASTEST:
      res = part::PartitionFastest(r, energy::BoltzEnergyModel::Create(em_));
      break;
    default: bug();
    }
    break;
  default: bug();
  }

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
  return {energy::EnergyModel::FromArgParse(args), CtxCfg::FromArgParse(args)};
}

}  // namespace mrna::ctx
