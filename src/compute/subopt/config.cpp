// Copyright 2022 E.
#include "compute/subopt/config.h"

#include "compute/energy/config.h"

namespace mrna::subopt {

void RegisterOpts(ArgParse* args) {
  energy::RegisterOpts(args);
  args->RegisterOpt(OPT_SUBOPT_DELTA);
  args->RegisterOpt(OPT_SUBOPT_MAX);
  args->RegisterOpt(OPT_SUBOPT_SORTED);
}

SuboptCfg SuboptCfg::FromArgParse(const ArgParse& args) {
  SuboptCfg cfg;
  if (args.Has(OPT_SUBOPT_DELTA)) cfg.delta = args.Get<Energy>(OPT_SUBOPT_DELTA);
  if (args.Has(OPT_SUBOPT_MAX)) cfg.strucs = args.Get<int>(OPT_SUBOPT_MAX);
  cfg.sorted = args.Has(OPT_SUBOPT_SORTED);
  return cfg;
}

}  // namespace mrna::subopt
