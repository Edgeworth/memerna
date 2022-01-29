// Copyright 2021 E.
#ifndef CONTEXT_CONFIG_H_
#define CONTEXT_CONFIG_H_

#include <string>

#include "util/argparse.h"

namespace mrna {

inline const Opt OPT_DP_ALG = Opt()
                                  .LongName("dp-alg")
                                  .Default("2")
                                  .Choice({"0", "1", "2", "3", "brute"})
                                  .Help("which algorithm for mfe folding");
inline const Opt OPT_SUBOPT_ALG = Opt()
                                      .LongName("subopt-alg")
                                      .Default("1")
                                      .Choice({"0", "1", "brute"})
                                      .Help("which algorithm for suboptimal folding");
inline const Opt OPT_PART_ALG = Opt()
                                    .LongName("part-alg")
                                    .Default("1")
                                    .Choice({"0", "1", "brute"})
                                    .Help("which algorithm for the partition function");

void RegisterOpts(ArgParse* args);

struct CtxCfg {
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

  CtxCfg(TableAlg table_alg_ = TableAlg::TWO, SuboptimalAlg suboptimal_alg_ = SuboptimalAlg::ONE,
      PartitionAlg partition_alg_ = PartitionAlg::ONE)
      : table_alg(table_alg_), suboptimal_alg(suboptimal_alg_), partition_alg(partition_alg_) {}

  TableAlg table_alg;
  SuboptimalAlg suboptimal_alg;
  PartitionAlg partition_alg;

  static CtxCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna

#endif  // CONTEXT_CONFIG_H_
