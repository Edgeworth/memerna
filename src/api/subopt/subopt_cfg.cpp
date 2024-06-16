// Copyright 2022 Eliot Courtney.
#include "api/subopt/subopt_cfg.h"

#include "api/ctx/backend_cfg.h"

namespace mrna::subopt {

void RegisterOpts(ArgParse* args) {
  RegisterOptsBackendCfg(args);
  args->RegisterOpt(OPT_SUBOPT_DELTA);
  args->RegisterOpt(OPT_SUBOPT_MAX);
  args->RegisterOpt(OPT_SUBOPT_TIME_SECS);
  args->RegisterOpt(OPT_SUBOPT_SORTED);
}

SuboptCfg SuboptCfg::FromArgParse(const ArgParse& args) {
  SuboptCfg cfg;
  args.MaybeSet(OPT_SUBOPT_DELTA, &cfg.delta);
  args.MaybeSet(OPT_SUBOPT_MAX, &cfg.strucs);
  args.MaybeSet(OPT_SUBOPT_TIME_SECS, &cfg.time_secs);
  args.MaybeSet(OPT_SUBOPT_SORTED, &cfg.sorted);
  return cfg;
}

}  // namespace mrna::subopt
