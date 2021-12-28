// Copyright 2021 E.
// Copyright 2021 E.
#include "model/config.h"

namespace mrna {

ModelCfg ModelCfgFromArgParse(const ArgParse& args) {
  ModelCfg cfg;
  const auto dp_alg = args.GetOption("dp-alg");
  if (dp_alg == "0") {
    cfg.table_alg = ModelCfg::TableAlg::ZERO;
  } else if (dp_alg == "1") {
    cfg.table_alg = ModelCfg::TableAlg::ONE;
  } else if (dp_alg == "2") {
    cfg.table_alg = ModelCfg::TableAlg::TWO;
  } else if (dp_alg == "3") {
    cfg.table_alg = ModelCfg::TableAlg::THREE;
  } else if (dp_alg == "brute") {
    cfg.table_alg = ModelCfg::TableAlg::BRUTE;
  } else {
    error("unknown fold option");
  }
  const auto subopt_alg = args.GetOption("subopt-alg");
  if (subopt_alg == "0") {
    cfg.suboptimal_alg = ModelCfg::SuboptimalAlg::ZERO;
  } else if (subopt_alg == "1") {
    cfg.suboptimal_alg = ModelCfg::SuboptimalAlg::ONE;
  } else if (subopt_alg == "brute") {
    cfg.suboptimal_alg = ModelCfg::SuboptimalAlg::BRUTE;
  } else {
    error("unknown suboptimal option");
  }
  const auto part_alg = args.GetOption("part-alg");
  if (part_alg == "0") {
    cfg.partition_alg = ModelCfg::PartitionAlg::ZERO;
  } else if (part_alg == "1") {
    cfg.partition_alg = ModelCfg::PartitionAlg::ONE;
  } else if (part_alg == "brute") {
    cfg.partition_alg = ModelCfg::PartitionAlg::BRUTE;
  } else {
    error("unknown partition option");
  }
  return cfg;
}

}  // namespace mrna
