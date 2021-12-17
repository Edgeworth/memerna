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

constexpr context_opt_t::TableAlg context_opt_t::TABLE_ALGS[];
constexpr context_opt_t::SuboptimalAlg context_opt_t::SUBOPTIMAL_ALGS[];
constexpr context_opt_t::PartitionAlg context_opt_t::COMPUTE_PARTITION_ALGS[];

context_opt_t ContextOptionsFromArgParse(const ArgParse& args) {
  context_opt_t options;
  const auto dp_alg = args.GetOption("dp-alg");
  if (dp_alg == "0") {
    options.table_alg = context_opt_t::TableAlg::ZERO;
  } else if (dp_alg == "1") {
    options.table_alg = context_opt_t::TableAlg::ONE;
  } else if (dp_alg == "2") {
    options.table_alg = context_opt_t::TableAlg::TWO;
  } else if (dp_alg == "3") {
    options.table_alg = context_opt_t::TableAlg::THREE;
  } else if (dp_alg == "brute") {
    options.table_alg = context_opt_t::TableAlg::BRUTE;
  } else {
    verify(false, "unknown fold option");
  }
  const auto subopt_alg = args.GetOption("subopt-alg");
  if (subopt_alg == "0") {
    options.suboptimal_alg = context_opt_t::SuboptimalAlg::ZERO;
  } else if (subopt_alg == "1") {
    options.suboptimal_alg = context_opt_t::SuboptimalAlg::ONE;
  } else if (subopt_alg == "brute") {
    options.suboptimal_alg = context_opt_t::SuboptimalAlg::BRUTE;
  } else {
    verify(false, "unknown suboptimal option");
  }
  const auto part_alg = args.GetOption("part-alg");
  if (part_alg == "0") {
    options.partition_alg = context_opt_t::PartitionAlg::ZERO;
  } else if (part_alg == "1") {
    options.partition_alg = context_opt_t::PartitionAlg::ONE;
  } else if (part_alg == "brute") {
    options.partition_alg = context_opt_t::PartitionAlg::BRUTE;
  } else {
    verify(false, "unknown partition option");
  }
  return options;
}

void Context::ComputeTables() {
  mfe::SetMfeGlobalState(r, *em);
  switch (options.table_alg) {
  case context_opt_t::TableAlg::ZERO: mfe::internal::ComputeTables0(); break;
  case context_opt_t::TableAlg::ONE: mfe::internal::ComputeTables1(); break;
  case context_opt_t::TableAlg::TWO: mfe::internal::ComputeTables2(); break;
  case context_opt_t::TableAlg::THREE: mfe::internal::ComputeTables3(); break;
  default: verify(false, "bug");
  }
  mfe::internal::ComputeExterior();
}

computed_t Context::Fold() {
  if (options.table_alg == context_opt_t::TableAlg::BRUTE) return mfe::MfeBruteForce(r, *em);

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
  if (options.suboptimal_alg == context_opt_t::SuboptimalAlg::BRUTE) {
    auto computeds = subopt::SuboptimalBruteForce(r, *em, subopt_num);
    for (const auto& computed : computeds) fn(computed);
    return static_cast<int>(computeds.size());
  }

  ComputeTables();
  switch (options.suboptimal_alg) {
  case context_opt_t::SuboptimalAlg::ZERO:
    return subopt::internal::Suboptimal0(subopt_delta, subopt_num).Run(fn);
  case context_opt_t::SuboptimalAlg::ONE:
    return subopt::internal::Suboptimal1(subopt_delta, subopt_num).Run(fn, sorted);
  default:
    verify(
        false, "bug - no such suboptimal algorithm %d", static_cast<int>(options.suboptimal_alg));
  }
}

partition::partition_t Context::Partition() {
  partition::SetPartitionGlobalState(r, *em);
  switch (options.partition_alg) {
  case context_opt_t::PartitionAlg::ZERO: partition::internal::Partition0(); break;
  case context_opt_t::PartitionAlg::ONE: partition::internal::Partition1(); break;
  case context_opt_t::PartitionAlg::BRUTE: return partition::PartitionBruteForce(r, *em).first;
  }
  const auto& gpt = partition::internal::gpt;
  const int size = static_cast<int>(r.size());
  array3d_t<penergy_t, 1> p((std::size_t(size)));
  for (int i = 0; i < size; ++i)  // TODO optimise this?
    for (int j = 0; j < size; ++j) p[i][j][0] = gpt[i][j][partition::PT_P];
  return {std::move(p), partition::internal::gptext[0][partition::PTEXT_R]};
}

}  // namespace mrna
