// Copyright 2021 E.
#ifndef MODEL_CONFIG_H_
#define MODEL_CONFIG_H_

#include <map>
#include <string>

#include "util/argparse.h"

namespace mrna {

inline const std::map<std::string, Opt> MODEL_OPTS = {
    {"dp-alg", Opt("which algorithm for mfe folding").Arg("2", {"0", "1", "2", "3", "brute"})},
    {"subopt-alg", Opt("which algorithm for suboptimal folding").Arg("1", {"0", "1", "brute"})},
    {"part-alg", Opt("which algorithm for the partition function").Arg("1", {"0", "1", "brute"})}};

struct ModelCfg {
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
  inline static constexpr PartitionAlg PARTITION_ALGS[] = {PartitionAlg::ZERO, PartitionAlg::ONE};

  ModelCfg(TableAlg table_alg_ = TableAlg::TWO, SuboptimalAlg suboptimal_alg_ = SuboptimalAlg::ONE,
      PartitionAlg partition_alg_ = PartitionAlg::ONE)
      : table_alg(table_alg_), suboptimal_alg(suboptimal_alg_), partition_alg(partition_alg_) {}

  TableAlg table_alg;
  SuboptimalAlg suboptimal_alg;
  PartitionAlg partition_alg;

  static ModelCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna

#endif  // MODEL_CONFIG_H_
