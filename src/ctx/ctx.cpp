// Copyright 2016 E.
#include "ctx/ctx.h"

#include <cassert>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include "compute/boltz_dp.h"
#include "compute/brute/alg.h"
#include "compute/dp.h"
#include "compute/energy/boltzmann_model.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt0.h"
#include "compute/subopt/subopt1.h"
#include "compute/traceback/traceback.h"
#include "model/primary.h"
#include "util/array.h"
#include "util/error.h"

namespace mrna::ctx {

energy::EnergyResult Ctx::Efn(
    Primary r, Secondary s, const Ctds* given_ctd, bool build_structure) const {
  return em_.TotalEnergy(r, s, given_ctd, build_structure);
}

DpArray Ctx::ComputeTables(const Primary& r) const {
  switch (cfg_.dp_alg) {
  case CtxCfg::DpAlg::ZERO: return mfe::ComputeTables0(r, em_);
  case CtxCfg::DpAlg::ONE: return mfe::ComputeTables1(r, em_);
  case CtxCfg::DpAlg::TWO: return mfe::ComputeTables2(r, em_);
  case CtxCfg::DpAlg::THREE: return mfe::ComputeTables3(r, em_);
  default: bug();
  }
}

ctx::FoldResult Ctx::Fold(Primary r) const {
  if (cfg_.dp_alg == CtxCfg::DpAlg::BRUTE) {
    auto subopt = brute::MfeBruteForce(std::move(r), em_);
    return {.mfe = {.dp{}, .ext{}, .energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  auto dp = ComputeTables(r);
  auto ext = mfe::ComputeExterior(r, em_, dp);
  auto tb = tb::Traceback(r, em_, dp, ext);
  auto energy = ext[0][EXT];
  return ctx::FoldResult{
      .mfe = {.dp = std::move(dp), .ext = std::move(ext), .energy = energy},
      .tb = std::move(tb),
  };
}

std::vector<subopt::SuboptResult> Ctx::SuboptimalIntoVector(
    Primary r, subopt::SuboptCfg cfg) const {
  std::vector<subopt::SuboptResult> subopts;
  [[maybe_unused]] int strucs = Suboptimal(
      std::move(r), [&subopts](const subopt::SuboptResult& subopt) { subopts.push_back(subopt); },
      cfg);
  assert(strucs == static_cast<int>(subopts.size()));
  return subopts;
}

int Ctx::Suboptimal(Primary r, subopt::SuboptCallback fn, subopt::SuboptCfg cfg) const {
  if (cfg_.subopt_alg == CtxCfg::SuboptAlg::BRUTE) {
    // TODO: handle cases other than max structures.
    auto subopts = brute::SuboptimalBruteForce(std::move(r), em_, std::move(cfg));
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  auto dp = ComputeTables(r);
  auto ext = mfe::ComputeExterior(r, em_, dp);
  switch (cfg_.subopt_alg) {
  case CtxCfg::SuboptAlg::ZERO:
    return subopt::Suboptimal0(std::move(r), em_, std::move(dp), std::move(ext), cfg).Run(fn);
  case CtxCfg::SuboptAlg::ONE:
    return subopt::Suboptimal1(std::move(r), em_, std::move(dp), std::move(ext), cfg).Run(fn);
  default:
    verify(false, "bug - no such suboptimal algorithm %d", static_cast<int>(cfg_.subopt_alg));
  }
}

part::PartResult Ctx::Partition(Primary r) const {
  std::tuple<BoltzDpArray, BoltzExtArray> res;
  switch (cfg_.part_alg) {
  case CtxCfg::PartAlg::ZERO: res = part::Partition0(r, em_); break;
  case CtxCfg::PartAlg::ONE: res = part::Partition1(r, energy::BoltzEnergyModel(em_)); break;
  case CtxCfg::PartAlg::BRUTE: return brute::PartitionBruteForce(std::move(r), em_);
  }
  const int N = static_cast<int>(r.size());
  auto [dp, ext] = std::move(res);
  BoltzSums p(N, 0);
  for (int i = 0; i < N; ++i)  // TODO optimise this?
    for (int j = 0; j < N; ++j) p[i][j] = dp[i][j][PT_P];

  part::Part part{std::move(p), ext[0][PTEXT_R]};
  auto prob = part.Prob();
  return part::PartResult{
      .dp = std::move(dp), .ext = std::move(ext), .part = std::move(part), .prob = std::move(prob)};
}

Ctx Ctx::FromArgParse(const ArgParse& args) {
  return Ctx(energy::EnergyModel::FromArgParse(args), CtxCfg::FromArgParse(args));
}

}  // namespace mrna::ctx
