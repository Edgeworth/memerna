// Copyright 2021 Eliot Courtney.
#include "context/config.h"

#include "compute/energy/energy.h"
#include "util/error.h"

namespace mrna {

void RegisterOpts(ArgParse* args) {
  energy::RegisterOpts(args);
  args->RegisterOpt(OPT_DP_ALG);
  args->RegisterOpt(OPT_SUBOPT_ALG);
  args->RegisterOpt(OPT_PART_ALG);
}

CtxCfg CtxCfg::FromArgParse(const ArgParse& args) {
  CtxCfg cfg;
  const auto dp_alg = args.Get(OPT_DP_ALG);
  if (dp_alg == "0") {
    cfg.table_alg = CtxCfg::TableAlg::ZERO;
  } else if (dp_alg == "1") {
    cfg.table_alg = CtxCfg::TableAlg::ONE;
  } else if (dp_alg == "2") {
    cfg.table_alg = CtxCfg::TableAlg::TWO;
  } else if (dp_alg == "3") {
    cfg.table_alg = CtxCfg::TableAlg::THREE;
  } else if (dp_alg == "brute") {
    cfg.table_alg = CtxCfg::TableAlg::BRUTE;
  } else {
    error("unknown fold option");
  }
  const auto subopt_alg = args.Get(OPT_SUBOPT_ALG);
  if (subopt_alg == "0") {
    cfg.suboptimal_alg = CtxCfg::SuboptimalAlg::ZERO;
  } else if (subopt_alg == "1") {
    cfg.suboptimal_alg = CtxCfg::SuboptimalAlg::ONE;
  } else if (subopt_alg == "brute") {
    cfg.suboptimal_alg = CtxCfg::SuboptimalAlg::BRUTE;
  } else {
    error("unknown suboptimal option");
  }
  const auto part_alg = args.Get(OPT_PART_ALG);
  if (part_alg == "0") {
    cfg.partition_alg = CtxCfg::PartitionAlg::ZERO;
  } else if (part_alg == "1") {
    cfg.partition_alg = CtxCfg::PartitionAlg::ONE;
  } else if (part_alg == "brute") {
    cfg.partition_alg = CtxCfg::PartitionAlg::BRUTE;
  } else {
    error("unknown partition option");
  }
  return cfg;
}

}  // namespace mrna
