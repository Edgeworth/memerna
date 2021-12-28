// Copyright 2021 Eliot Courtney.
#ifndef MODEL_CONFIG_H_
#define MODEL_CONFIG_H_

#include "util/argparse.h"

namespace mrna {

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

  inline static constexpr TableAlg TABLE_ALGS[] = {
      TableAlg::ZERO, TableAlg::ONE, TableAlg::TWO, TableAlg::THREE};
  inline static constexpr SuboptimalAlg SUBOPTIMAL_ALGS[] = {
      SuboptimalAlg::ZERO, SuboptimalAlg::ONE};
  inline static constexpr PartitionAlg COMPUTE_PARTITION_ALGS[] = {
      PartitionAlg::ZERO, PartitionAlg::ONE};

  context_opt_t(TableAlg table_alg_ = TableAlg::ZERO,
      SuboptimalAlg suboptimal_alg_ = SuboptimalAlg::ZERO,
      PartitionAlg partition_alg_ = PartitionAlg::ZERO)
      : table_alg(table_alg_), suboptimal_alg(suboptimal_alg_), partition_alg(partition_alg_) {}

  TableAlg table_alg;
  SuboptimalAlg suboptimal_alg;
  PartitionAlg partition_alg;
};

inline const std::map<std::string, opt_t> CONTEXT_OPTIONS = {
    {"dp-alg", opt_t("which algorithm for mfe folding").Arg("2", {"0", "1", "2", "3", "brute"})},
    {"subopt-alg", opt_t("which algorithm for suboptimal folding").Arg("1", {"0", "1", "brute"})},
    {"part-alg",
        opt_t("which algorithm for the partition function").Arg("1", {"0", "1", "brute"})}};

context_opt_t ContextOptionsFromArgParse(const ArgParse& args);

}  // namespace mrna

#endif  // MODEL_CONFIG_H_
