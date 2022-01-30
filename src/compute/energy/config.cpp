// Copyright 2022 E.
#include "compute/energy/config.h"

#include "util/error.h"

namespace mrna::energy {

void RegisterOpts(ArgParse* args) {
  args->RegisterOpt(OPT_SEED);
  args->RegisterOpt(OPT_MEMERNA_DATA);
  args->RegisterOpt(OPT_LONELY_PAIRS);
  args->RegisterOpt(OPT_CTD);
}

EnergyCfg EnergyCfg::FromArgParse(const ArgParse& args) {
  EnergyCfg cfg;
  cfg.lonely_pairs = args.Has(OPT_LONELY_PAIRS);
  cfg.ctd = args.Get<Ctd>(OPT_CTD);
  return cfg;
}

std::istream& operator>>(std::istream& str, EnergyCfg::Ctd& o) {
  std::string s;
  str >> s;
  if (s == "none")
    o = EnergyCfg::Ctd::NONE;
  else if (s == "no_coax")
    o = EnergyCfg::Ctd::NO_COAX;
  else if (s == "all")
    o = EnergyCfg::Ctd::ALL;
  else
    error("Invalid CTD option %s", s.c_str());
  return str;
}

}  // namespace mrna::energy
