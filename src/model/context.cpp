// Copyright 2016 Eliot Courtney.
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
  mfe::SetMfeGlobalState(r, *em);
  switch (cfg.table_alg) {
  case ModelCfg::TableAlg::ZERO: mfe::internal::ComputeTables0(); break;
  case ModelCfg::TableAlg::ONE: mfe::internal::ComputeTables1(); break;
  case ModelCfg::TableAlg::TWO: mfe::internal::ComputeTables2(); break;
  case ModelCfg::TableAlg::THREE: mfe::internal::ComputeTables3(); break;
  default: bug();
  }
  mfe::internal::ComputeExterior();
}

computed_t Context::Fold() {
  if (cfg.table_alg == ModelCfg::TableAlg::BRUTE) return mfe::MfeBruteForce(r, *em);

  ComputeTables();
  traceback::Traceback();
  return {{gr, mfe::internal::gp}, mfe::internal::gctd, mfe::internal::genergy};
}

std::vector<computed_t> Context::SuboptimalIntoVector(
    bool sorted, energy_t subopt_delta, int subopt_num) {
  std::vector<computed_t> computeds;
  [[maybe_unused]] int num_structures =
      Suboptimal([&computeds](const computed_t& c) { computeds.push_back(c); }, sorted,
          subopt_delta, subopt_num);
  assert(num_structures == static_cast<int>(computeds.size()));
  return computeds;
}

int Context::Suboptimal(
    subopt::SuboptimalCallback fn, bool sorted, energy_t subopt_delta, int subopt_num) {
  if (cfg.suboptimal_alg == ModelCfg::SuboptimalAlg::BRUTE) {
    auto computeds = subopt::SuboptimalBruteForce(r, *em, subopt_num);
    for (const auto& computed : computeds) fn(computed);
    return static_cast<int>(computeds.size());
  }

  ComputeTables();
  switch (cfg.suboptimal_alg) {
  case ModelCfg::SuboptimalAlg::ZERO:
    return subopt::internal::Suboptimal0(subopt_delta, subopt_num).Run(fn);
  case ModelCfg::SuboptimalAlg::ONE:
    return subopt::internal::Suboptimal1(subopt_delta, subopt_num).Run(fn, sorted);
  default:
    verify(
        false, "bug - no such suboptimal algorithm %d", static_cast<int>(cfg.suboptimal_alg));
  }
}

partition::partition_t Context::Partition() {
  partition::SetPartitionGlobalState(r, *em);
  switch (cfg.partition_alg) {
  case ModelCfg::PartitionAlg::ZERO: partition::internal::Partition0(); break;
  case ModelCfg::PartitionAlg::ONE: partition::internal::Partition1(); break;
  case ModelCfg::PartitionAlg::BRUTE: return partition::PartitionBruteForce(r, *em).first;
  }
  const auto& gpt = partition::internal::gpt;
  const int size = static_cast<int>(r.size());
  array3d_t<penergy_t, 1> p((std::size_t(size)));
  for (int i = 0; i < size; ++i)  // TODO optimise this?
    for (int j = 0; j < size; ++j) p[i][j][0] = gpt[i][j][partition::PT_P];
  return {std::move(p), partition::internal::gptext[0][partition::PTEXT_R]};
}

}  // namespace mrna
