// Copyright 2016 Eliot Courtney.
#include "model/context.h"

#include <utility>
#include <vector>

#include "compute/dp.h"
#include "compute/mfe/brute.h"
#include "compute/partition/brute.h"
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

// TODO: Come up with better types than this tuple.
std::tuple<Computed, DpArray> Context::Fold() {
  if (cfg_.table_alg == ModelCfg::TableAlg::BRUTE) return {mfe::MfeBruteForce(r_, em_), DpArray()};

  auto dp = ComputeTables();
  auto ext = mfe::ComputeExterior(r_, em_, dp);
  auto [s, ctd] = traceback::Traceback(r_, em_, dp, ext);
  return {Computed{r_, std::move(s), std::move(ctd), ext[0][EXT]}, std::move(dp)};
}

std::vector<Computed> Context::SuboptimalIntoVector(
    bool sorted, Energy subopt_delta, int subopt_num) {
  std::vector<Computed> computeds;
  [[maybe_unused]] int num_structures =
      Suboptimal([&computeds](const Computed& c) { computeds.push_back(c); }, sorted, subopt_delta,
          subopt_num);
  assert(num_structures == static_cast<int>(computeds.size()));
  return computeds;
}

int Context::Suboptimal(
    subopt::SuboptimalCallback fn, bool sorted, Energy subopt_delta, int subopt_num) {
  if (cfg_.suboptimal_alg == ModelCfg::SuboptimalAlg::BRUTE) {
    auto computeds = subopt::SuboptimalBruteForce(r_, em_, subopt_num);
    for (const auto& computed : computeds) fn(computed);
    return static_cast<int>(computeds.size());
  }

  auto dp = ComputeTables();
  auto ext = mfe::ComputeExterior(r_, em_, dp);
  switch (cfg_.suboptimal_alg) {
  case ModelCfg::SuboptimalAlg::ZERO:
    return subopt::Suboptimal0(r_, em_, std::move(dp), std::move(ext), subopt_delta, subopt_num)
        .Run(fn);
  case ModelCfg::SuboptimalAlg::ONE:
    return subopt::Suboptimal1(r_, em_, std::move(dp), std::move(ext), subopt_delta, subopt_num)
        .Run(fn, sorted);
  default:
    verify(false, "bug - no such suboptimal algorithm %d", static_cast<int>(cfg_.suboptimal_alg));
  }
}

partition::Partition Context::Partition() {
  std::tuple<BoltzDpArray, BoltzExtArray> res;
  switch (cfg_.partition_alg) {
  case ModelCfg::PartitionAlg::ZERO: res = partition::Partition0(r_, em_); break;
  case ModelCfg::PartitionAlg::ONE:
    res = partition::Partition1(r_, energy::BoltzEnergyModel(em_));
    break;
  case ModelCfg::PartitionAlg::BRUTE: return partition::PartitionBruteForce(r_, em_).first;
  }
  const int N = static_cast<int>(r_.size());
  auto [dp, ext] = std::move(res);
  Array3D<BoltzEnergy, 1> p(N);
  for (int i = 0; i < N; ++i)  // TODO optimise this?
    for (int j = 0; j < N; ++j) p[i][j][0] = dp[i][j][PT_P];
  return {std::move(p), ext[0][PTEXT_R]};
}

}  // namespace mrna
