// Copyright 2016 Eliot Courtney.
#ifndef MEMERNA_CONTEXT_H_
#define MEMERNA_CONTEXT_H_

#include "argparse.h"
#include "common.h"
#include "energy/energy.h"
#include "fold/fold.h"
#include "partition/partition.h"

namespace memerna {

struct context_opt_t {
  enum class TableAlg {
    ZERO,
    ONE,
    TWO,
    THREE,
    BRUTE  // Not included in the normal table algs since exponential.
  };

  enum class SuboptimalAlg { ZERO, ONE, BRUTE };

  enum class PartitionAlg { ZERO, ONE, BRUTE };

  static constexpr TableAlg TABLE_ALGS[] = {
      TableAlg::ZERO, TableAlg::ONE, TableAlg::TWO, TableAlg::THREE};
  static constexpr SuboptimalAlg SUBOPTIMAL_ALGS[] = {SuboptimalAlg::ZERO, SuboptimalAlg::ONE};
  static constexpr PartitionAlg PARTITION_ALGS[] = {PartitionAlg::ZERO, PartitionAlg::ONE};

  context_opt_t(TableAlg table_alg_ = TableAlg::ZERO,
      SuboptimalAlg suboptimal_alg_ = SuboptimalAlg::ZERO,
      PartitionAlg partition_alg_ = PartitionAlg::ZERO)
      : table_alg(table_alg_), suboptimal_alg(suboptimal_alg_), partition_alg(partition_alg_) {}

  TableAlg table_alg;
  SuboptimalAlg suboptimal_alg;
  PartitionAlg partition_alg;
};

class Context {
public:
  Context(const primary_t& r_, const energy::EnergyModelPtr em_) : r(r_), em(em_), options() {
    verify_expr(r.size() > 0u, "cannot process zero length RNA");
  };
  Context(const primary_t& r_, const energy::EnergyModelPtr em_, context_opt_t options_)
      : r(r_), em(em_), options(options_) {
    verify_expr(r.size() > 0u, "cannot process zero length RNA");
  }
  Context(const Context& o) = default;

  Context() = delete;
  Context& operator=(const Context&) = delete;
  Context(Context&& o) = delete;
  Context& operator=(Context&&) = delete;

  computed_t Fold();
  std::vector<computed_t> SuboptimalIntoVector(
      bool sorted, energy_t subopt_delta = -1, int subopt_num = -1);
  int Suboptimal(
      fold::SuboptimalCallback fn, bool sorted, energy_t subopt_delta = -1, int subopt_num = -1);
  partition::partition_t Partition();

private:
  const primary_t r;
  const energy::EnergyModelPtr em;
  const context_opt_t options;

  void ComputeTables();
};

const std::map<std::string, opt_t> CONTEXT_OPTIONS = {
    {"dp-alg", opt_t("which algorithm for mfe folding").Arg("2", {"0", "1", "2", "3", "brute"})},
    {"subopt-alg", opt_t("which algorithm for suboptimal folding").Arg("1", {"0", "1", "brute"})},
    {"part-alg",
        opt_t("which algorithm for the partition function").Arg("1", {"0", "1", "brute"})}};

context_opt_t ContextOptionsFromArgParse(const ArgParse& argparse);

}  // namespace memerna

#endif  // MEMERNA_CONTEXT_H_
