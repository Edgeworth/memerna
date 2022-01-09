// Copyright 2016 E.
#include "model/context.h"

#include <utility>
#include <vector>

#include "compute/mfe/brute.h"
#include "compute/partition/brute.h"
#include "compute/partition/globals.h"
#include "compute/subopt/brute.h"
#include "compute/subopt/subopt0.h"
#include "compute/subopt/subopt1.h"
#include "compute/traceback/traceback.h"

namespace mrna {

void Context::ComputeTables() {
  mfe::SetMfeGlobalState(r_);
  switch (cfg_.table_alg) {
  case ModelCfg::TableAlg::ZERO: mfe::internal::ComputeTables0(r_, em_); break;
  case ModelCfg::TableAlg::ONE: mfe::internal::ComputeTables1(r_, em_); break;
  case ModelCfg::TableAlg::TWO: mfe::internal::ComputeTables2(r_, em_); break;
  case ModelCfg::TableAlg::THREE: mfe::internal::ComputeTables3(r_, em_); break;
  default: bug();
  }
  mfe::internal::ComputeExterior(r_, em_);
}

Computed Context::Fold() {
  if (cfg_.table_alg == ModelCfg::TableAlg::BRUTE) return mfe::MfeBruteForce(r_, em_);

  ComputeTables();
  traceback::Traceback(r_, em_);
  return {{r_, mfe::internal::gp}, mfe::internal::gctd, mfe::internal::genergy};
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

  ComputeTables();
  switch (cfg_.suboptimal_alg) {
  case ModelCfg::SuboptimalAlg::ZERO:
    return subopt::internal::Suboptimal0(r_, em_, subopt_delta, subopt_num).Run(fn);
  case ModelCfg::SuboptimalAlg::ONE:
    return subopt::internal::Suboptimal1(r_, em_, subopt_delta, subopt_num).Run(fn, sorted);
  default:
    verify(false, "bug - no such suboptimal algorithm %d", static_cast<int>(cfg_.suboptimal_alg));
  }
}

partition::Partition Context::Partition() {
  partition::SetPartitionGlobalState(r_);
  switch (cfg_.partition_alg) {
  case ModelCfg::PartitionAlg::ZERO: partition::internal::Partition0(r_, em_); break;
  case ModelCfg::PartitionAlg::ONE:
    partition::internal::Partition1(r_, energy::BoltzEnergyModel(em_));
    break;
  case ModelCfg::PartitionAlg::BRUTE: return partition::PartitionBruteForce(r_, em_).first;
  }
  const auto& gpt = partition::internal::gpt;
  const int size = static_cast<int>(r_.size());
  Array3D<BoltzEnergy, 1> p((std::size_t(size)));
  for (int i = 0; i < size; ++i)  // TODO optimise this?
    for (int j = 0; j < size; ++j) p[i][j][0] = gpt[i][j][partition::PT_P];
  return {std::move(p), partition::internal::gptext[0][partition::PTEXT_R]};
}

}  // namespace mrna
