// Copyright 2022 Eliot Courtney.
#include "compute/brute/config.h"

#include "options.h"

namespace mrna::brute {

void RegisterOpts(ArgParse* args) {
  subopt::RegisterOpts(args);
  args->RegisterOpt(OPT_FOLD);
  args->RegisterOpt(OPT_SUBOPT);
  args->RegisterOpt(OPT_PART);
}

BruteCfg BruteCfg::FromArgParse(const ArgParse& args) {
  BruteCfg cfg;
  args.MaybeSet(OPT_FOLD, &cfg.mfe);
  args.MaybeSet(OPT_SUBOPT, &cfg.subopt);
  args.MaybeSet(OPT_PART, &cfg.part);
  cfg.subopt_cfg = subopt::SuboptCfg::FromArgParse(args);
  return cfg;
}

}  // namespace mrna::brute
