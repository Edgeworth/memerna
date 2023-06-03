// Copyright 2022 Eliot Courtney.
#include "api/subopt/subopt_cfg.h"

#include "api/energy/energy_cfg.h"
#include "api/energy/model.h"

namespace mrna::subopt {

void RegisterOpts(ArgParse* args) {
  erg::RegisterOptsEnergyModel(args);
  args->RegisterOpt(OPT_SUBOPT_DELTA);
  args->RegisterOpt(OPT_SUBOPT_MAX);
  args->RegisterOpt(OPT_SUBOPT_SORTED);
}

SuboptCfg SuboptCfg::FromArgParse(const ArgParse& args) {
  SuboptCfg cfg;
  args.MaybeSet(OPT_SUBOPT_DELTA, &cfg.delta);
  args.MaybeSet(OPT_SUBOPT_MAX, &cfg.strucs);
  args.MaybeSet(OPT_SUBOPT_SORTED, &cfg.sorted);
  return cfg;
}

}  // namespace mrna::subopt
