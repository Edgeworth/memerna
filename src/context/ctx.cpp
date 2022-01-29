// Copyright 2016 Eliot Courtney.
#include "context/ctx.h"

#include <cassert>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include "compute/boltz_dp.h"
#include "compute/dp.h"
#include "compute/energy/boltzmann_model.h"
#include "compute/mfe/brute.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/brute.h"
#include "compute/partition/partition.h"
#include "compute/subopt/brute.h"
#include "compute/subopt/subopt0.h"
#include "compute/subopt/subopt1.h"
#include "compute/traceback/traceback.h"
#include "model/primary.h"
#include "util/array.h"
#include "util/error.h"

namespace mrna {

DpArray Ctx::ComputeTables(const Primary& r) const {
  switch (cfg_.table_alg) {
  case CtxCfg::TableAlg::ZERO: return mfe::ComputeTables0(r, em_);
  case CtxCfg::TableAlg::ONE: return mfe::ComputeTables1(r, em_);
  case CtxCfg::TableAlg::TWO: return mfe::ComputeTables2(r, em_);
  case CtxCfg::TableAlg::THREE: return mfe::ComputeTables3(r, em_);
  default: bug();
  }
}

FoldResult Ctx::Fold(Primary r) const {
  if (cfg_.table_alg == CtxCfg::TableAlg::BRUTE) {
    auto subopt = mfe::MfeBruteForce(std::move(r), em_);
    return {
        .mfe = mfe::MfeResult{.dp{}, .ext{}, .energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  auto dp = ComputeTables(r);
  auto ext = mfe::ComputeExterior(r, em_, dp);
  auto tb = tb::Traceback(r, em_, dp, ext);
  auto energy = ext[0][EXT];
  return FoldResult{
      .mfe = mfe::MfeResult{.dp = std::move(dp), .ext = std::move(ext), .energy = energy},
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
  if (cfg_.suboptimal_alg == CtxCfg::SuboptimalAlg::BRUTE) {
    // TODO: handle cases other than max structures.
    auto subopts = subopt::SuboptimalBruteForce(std::move(r), em_, cfg.strucs);
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  auto dp = ComputeTables(r);
  auto ext = mfe::ComputeExterior(r, em_, dp);
  switch (cfg_.suboptimal_alg) {
  case CtxCfg::SuboptimalAlg::ZERO:
    return subopt::Suboptimal0(std::move(r), em_, std::move(dp), std::move(ext), cfg).Run(fn);
  case CtxCfg::SuboptimalAlg::ONE:
    return subopt::Suboptimal1(std::move(r), em_, std::move(dp), std::move(ext), cfg).Run(fn);
  default:
    verify(false, "bug - no such suboptimal algorithm %d", static_cast<int>(cfg_.suboptimal_alg));
  }
}

partition::PartitionResult Ctx::Partition(Primary r) const {
  std::tuple<BoltzDpArray, BoltzExtArray> res;
  switch (cfg_.partition_alg) {
  case CtxCfg::PartitionAlg::ZERO: res = partition::Partition0(r, em_); break;
  case CtxCfg::PartitionAlg::ONE:
    res = partition::Partition1(r, energy::BoltzEnergyModel(em_));
    break;
  case CtxCfg::PartitionAlg::BRUTE: return partition::PartitionBruteForce(std::move(r), em_);
  }
  const int N = static_cast<int>(r.size());
  auto [dp, ext] = std::move(res);
  BoltzSums p(N, 0);
  for (int i = 0; i < N; ++i)  // TODO optimise this?
    for (int j = 0; j < N; ++j) p[i][j] = dp[i][j][PT_P];

  partition::Partition part{std::move(p), ext[0][PTEXT_R]};
  auto prob = partition::ComputeBoltzProbs(part);
  return partition::PartitionResult{
      .dp = std::move(dp), .ext = std::move(ext), .p = std::move(part), .prob = std::move(prob)};
}

Ctx Ctx::FromArgParse(const ArgParse& args) {
  return Ctx(energy::EnergyModel::FromArgParse(args), CtxCfg::FromArgParse(args));
}

}  // namespace mrna
