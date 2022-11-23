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
#include "compute/partition/t04/partition.h"
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

DpArray Ctx::ComputeMfe(const Primary& r) const {
  return std::visit(
      [&](const energy::t04::ModelPtr& em) -> DpArray {
        switch (cfg_.dp_alg) {
        case CtxCfg::DpAlg::SLOWEST: return mfe::t04::MfeSlowest(r, em);
        case CtxCfg::DpAlg::SLOW: return mfe::t04::MfeSlow(r, em);
        case CtxCfg::DpAlg::FASTEST: return mfe::t04::MfeFastest(r, em);
        case CtxCfg::DpAlg::LYNGSO: return mfe::t04::MfeLyngso(r, em);
        default: bug();
        }
      },
      em_);
}

std::tuple<ExtArray, Energy> Ctx::ComputeMfeExterior(const Primary& r, const DpArray& dp) const {
  return std::visit(
      [&](const energy::t04::ModelPtr& em) -> std::tuple<ExtArray, Energy> {
        auto ext = mfe::t04::MfeExterior(r, em, dp);
        auto energy = ext[0][EXT];
        return {std::move(ext), energy};
      },
      em_);
}

tb::TracebackResult Ctx::ComputeTraceback(
    const Primary& r, const DpArray& dp, const ExtArray& ext) const {
  return std::visit(
      [&](const energy::t04::ModelPtr& em) -> tb::TracebackResult {
        return tb::t04::Traceback(r, em, dp, ext);
      },
      em_);
}

ctx::FoldResult Ctx::Fold(const Primary& r) const {
  if (cfg_.dp_alg == CtxCfg::DpAlg::BRUTE) {
    auto subopt = brute::MfeBrute(r, em_);
    return {.mfe = {.dp{}, .ext{}, .energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  auto dp = ComputeMfe(r);
  auto [ext, energy] = ComputeMfeExterior(r, dp);
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
    auto subopts = brute::SuboptBrute(r, em_, cfg);
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  auto dp = ComputeMfe(r);
  auto [ext, energy] = ComputeMfeExterior(r, dp);

  return std::visit(
      [&, ext = std::move(ext)](const energy::t04::ModelPtr& em) mutable -> int {
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
      em_);
}

part::PartResult Ctx::Partition(const Primary& r) const {
  if (cfg_.dp_alg == CtxCfg::DpAlg::BRUTE) {
    return brute::PartitionBrute(r, em_);
  }

  std::tuple<BoltzDpArray, BoltzExtArray> res;
  std::visit(
      [&](const energy::t04::ModelPtr& em) {
        switch (cfg_.part_alg) {
        case CtxCfg::PartAlg::SLOWEST: res = part::t04::PartitionSlowest(r, em); break;
        case CtxCfg::PartAlg::FASTEST:
          res = part::t04::PartitionFastest(r, energy::t04::BoltzModel::Create(em));
          break;
        default: bug();
        }
      },
      em_);

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
  return {energy::FromArgParse(args), CtxCfg::FromArgParse(args)};
}

}  // namespace mrna::ctx
