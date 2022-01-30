// Copyright 2022 E.
#include "compute/brute/config.h"

#include "options.h"

namespace mrna::brute {

void RegisterOpts(ArgParse* args) {
  subopt::RegisterOpts(args);
  args->RegisterOpt(OPT_MFE);
  args->RegisterOpt(OPT_SUBOPT);
  args->RegisterOpt(OPT_PART);
}

BruteCfg BruteCfg::FromArgParse(const ArgParse& args) {
  BruteCfg cfg;
  cfg.mfe = args.Has(OPT_MFE);
  cfg.subopt = args.Has(OPT_SUBOPT);
  cfg.part = args.Has(OPT_PART);
  cfg.subopt_cfg = subopt::SuboptCfg::FromArgParse(args);
  return cfg;
}

}  // namespace mrna::brute
