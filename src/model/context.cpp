// Copyright 2016 E.
#include "model/context.h"

#include <utility>
#include <vector>

#include "compute/dp.h"
#include "compute/mfe/brute.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/brute.h"
#include "compute/partition/partition.h"
#include "compute/subopt/brute.h"
#include "compute/subopt/subopt0.h"
#include "compute/subopt/subopt1.h"
#include "compute/traceback/traceback.h"

namespace mrna {

DpArray Context::ComputeTables() {
  switch (cfg_.table_alg) {
  case ModelCfg::TableAlg::ZERO: return mfe::ComputeTables0(r_, em_);
  case ModelCfg::TableAlg::ONE: return mfe::ComputeTables1(r_, em_);
  case ModelCfg::TableAlg::TWO: return mfe::ComputeTables2(r_, em_);
  case ModelCfg::TableAlg::THREE: return mfe::ComputeTables3(r_, em_);
  default: bug();
  }
}

FoldResult Context::Fold() {
  if (cfg_.table_alg == ModelCfg::TableAlg::BRUTE) {
    auto subopt = mfe::MfeBruteForce(Primary(r_), em_);
    return {.mfe = mfe::MfeResult{.energy = subopt.energy}, .tb = std::move(subopt.tb)};
  }

  auto dp = ComputeTables();
  auto ext = mfe::ComputeExterior(r_, em_, dp);
  auto tb = tb::Traceback(r_, em_, dp, ext);
  return FoldResult{
      .mfe = mfe::MfeResult{.dp = std::move(dp), .ext = std::move(ext), .energy = ext[0][EXT]},
      .tb = std::move(tb),
  };
}

std::vector<subopt::SuboptResult> Context::SuboptimalIntoVector(
    bool sorted, Energy subopt_delta, int subopt_num) {
  std::vector<subopt::SuboptResult> subopts;
  [[maybe_unused]] int num_structures =
      Suboptimal([&subopts](const subopt::SuboptResult& subopt) { subopts.push_back(subopt); },
          sorted, subopt_delta, subopt_num);
  assert(num_structures == static_cast<int>(subopts.size()));
  return subopts;
}

int Context::Suboptimal(
    subopt::SuboptCallback fn, bool sorted, Energy subopt_delta, int subopt_num) {
  if (cfg_.suboptimal_alg == ModelCfg::SuboptimalAlg::BRUTE) {
    auto subopts = subopt::SuboptimalBruteForce(Primary(r_), em_, subopt_num);
    for (const auto& subopt : subopts) fn(subopt);
    return static_cast<int>(subopts.size());
  }

  auto dp = ComputeTables();
  auto ext = mfe::ComputeExterior(r_, em_, dp);
  switch (cfg_.suboptimal_alg) {
  case ModelCfg::SuboptimalAlg::ZERO:
    return subopt::Suboptimal0(
        Primary(r_), em_, std::move(dp), std::move(ext), subopt_delta, subopt_num)
        .Run(fn);
  case ModelCfg::SuboptimalAlg::ONE:
    return subopt::Suboptimal1(
        Primary(r_), em_, std::move(dp), std::move(ext), subopt_delta, subopt_num)
        .Run(fn, sorted);
  default:
    verify(false, "bug - no such suboptimal algorithm %d", static_cast<int>(cfg_.suboptimal_alg));
  }
}

partition::PartitionResult Context::Partition() {
  std::tuple<BoltzDpArray, BoltzExtArray> res;
  switch (cfg_.partition_alg) {
  case ModelCfg::PartitionAlg::ZERO: res = partition::Partition0(r_, em_); break;
  case ModelCfg::PartitionAlg::ONE:
    res = partition::Partition1(r_, energy::BoltzEnergyModel(em_));
    break;
  case ModelCfg::PartitionAlg::BRUTE: return partition::PartitionBruteForce(Primary(r_), em_);
  }
  const int N = static_cast<int>(r_.size());
  auto [dp, ext] = std::move(res);
  Array3D<BoltzEnergy, 1> p(N);
  for (int i = 0; i < N; ++i)  // TODO optimise this?
    for (int j = 0; j < N; ++j) p[i][j][0] = dp[i][j][PT_P];

  partition::Partition part{std::move(p), ext[0][PTEXT_R]};
  auto prob = partition::ComputeProbabilities(part);
  return partition::PartitionResult{.p = std::move(part), .prob = std::move(prob)};
}

}  // namespace mrna
