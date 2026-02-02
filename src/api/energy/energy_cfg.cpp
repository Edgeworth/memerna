// Copyright 2022 Eliot Courtney.
#include "api/energy/energy_cfg.h"

#include <ostream>

namespace mrna::erg {

void RegisterOptsEnergyCfg(ArgParse* args) {
  args->RegisterOpt(OPT_LONELY_PAIRS);
  args->RegisterOpt(OPT_BULGE_STATES);
  args->RegisterOpt(OPT_CTD);
}

EnergyCfg EnergyCfg::FromArgParse(const ArgParse& args) {
  EnergyCfg cfg;
  args.MaybeSet(OPT_LONELY_PAIRS, &cfg.lonely_pairs);
  args.MaybeSet(OPT_BULGE_STATES, &cfg.bulge_states);
  args.MaybeSet(OPT_CTD, &cfg.ctd);
  return cfg;
}

std::ostream& operator<<(std::ostream& str, const EnergyCfg& o) {
  return str << "EnergyCfg{"
             << "lonely_pairs=" << Conv(o.lonely_pairs) << ", ctd=" << Conv(o.ctd)
             << ", bulge_states=" << o.bulge_states << "}";
}

}  // namespace mrna::erg
