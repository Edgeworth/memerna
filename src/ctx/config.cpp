// Copyright 2021 E.
#include "ctx/config.h"

#include "compute/energy/config.h"
#include "compute/subopt/config.h"
#include "util/error.h"

namespace mrna::ctx {

void RegisterOpts(ArgParse* args) {
  energy::RegisterOpts(args);
  subopt::RegisterOpts(args);
  args->RegisterOpt(OPT_DP_ALG);
  args->RegisterOpt(OPT_SUBOPT_ALG);
  args->RegisterOpt(OPT_PART_ALG);
}

CtxCfg CtxCfg::FromArgParse(const ArgParse& args) {
  CtxCfg cfg;
  const auto dp_alg = args.Get(OPT_DP_ALG);
  if (dp_alg == "0") {
    cfg.dp_alg = CtxCfg::DpAlg::ZERO;
  } else if (dp_alg == "1") {
    cfg.dp_alg = CtxCfg::DpAlg::ONE;
  } else if (dp_alg == "2") {
    cfg.dp_alg = CtxCfg::DpAlg::TWO;
  } else if (dp_alg == "3") {
    cfg.dp_alg = CtxCfg::DpAlg::THREE;
  } else if (dp_alg == "brute") {
    cfg.dp_alg = CtxCfg::DpAlg::BRUTE;
  } else {
    error("unknown fold option");
  }
  const auto subopt_alg = args.Get(OPT_SUBOPT_ALG);
  if (subopt_alg == "0") {
    cfg.subopt_alg = CtxCfg::SuboptAlg::ZERO;
  } else if (subopt_alg == "1") {
    cfg.subopt_alg = CtxCfg::SuboptAlg::ONE;
  } else if (subopt_alg == "brute") {
    cfg.subopt_alg = CtxCfg::SuboptAlg::BRUTE;
  } else {
    error("unknown suboptimal option");
  }
  const auto part_alg = args.Get(OPT_PART_ALG);
  if (part_alg == "0") {
    cfg.part_alg = CtxCfg::PartAlg::ZERO;
  } else if (part_alg == "1") {
    cfg.part_alg = CtxCfg::PartAlg::ONE;
  } else if (part_alg == "brute") {
    cfg.part_alg = CtxCfg::PartAlg::BRUTE;
  } else {
    error("unknown partition option");
  }
  return cfg;
}

}  // namespace mrna::ctx
