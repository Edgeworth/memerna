// Copyright 2022 E.
#include "api/brute/brute_cfg.h"

#include "api/options.h"

namespace mrna::brute {

void RegisterOpts(ArgParse* args) {
  subopt::RegisterOpts(args);
  args->RegisterOpt(OPT_FOLD);
  args->RegisterOpt(OPT_SUBOPT);
  args->RegisterOpt(OPT_PFN);
}

BruteCfg BruteCfg::FromArgParse(const ArgParse& args) {
  BruteCfg cfg;
  args.MaybeSet(OPT_FOLD, &cfg.mfe);
  args.MaybeSet(OPT_SUBOPT, &cfg.subopt);
  args.MaybeSet(OPT_PFN, &cfg.pfn);
  cfg.subopt_cfg = subopt::SuboptCfg::FromArgParse(args);
  return cfg;
}

}  // namespace mrna::brute
