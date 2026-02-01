// Copyright 2022 Eliot Courtney.
#include "api/energy/energy_cfg.h"

#include <ostream>
#include <string>

#include "util/error.h"
#include "util/util.h"

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

void EnergyCfgSupport::VerifySupported(const std::string& name, const EnergyCfg& cfg) const {
  verify(Contains(lonely_pairs, cfg.lonely_pairs), "{} does not support lonely pairs option: {}",
      name, cfg.lonely_pairs);
  verify(Contains(bulge_states, cfg.bulge_states), "{} does not support bulge states option: {}",
      name, cfg.bulge_states);
  verify(Contains(ctd, cfg.ctd), "{} does not support CTD option: {}", name, cfg.ctd);
}

std::ostream& operator<<(std::ostream& str, const EnergyCfg& o) {
  return str << "EnergyCfg{"
             << "lonely_pairs=" << Conv(o.lonely_pairs) << ", ctd=" << Conv(o.ctd)
             << ", bulge_states=" << o.bulge_states << "}";
}

}  // namespace mrna::erg
