// Copyright 2022 E.
#include "model/config.h"

namespace mrna {

void RegisterOpts(ArgParse* args) {
  args->RegisterOpt(OPT_LONELY_PAIRS);
  args->RegisterOpt(OPT_CTD);
}

ModelCfg ModelCfg::FromArgParse(const ArgParse& args) {
  ModelCfg cfg;
  cfg.lonely_pairs = args.Has(OPT_LONELY_PAIRS);
  cfg.ctd = args.Get<Ctd>(OPT_CTD);
  return cfg;
}

std::istream& operator>>(std::istream& str, ModelCfg::Ctd& o) {
  std::string s;
  str >> s;
  if (s == "none")
    o = ModelCfg::Ctd::NONE;
  else if (s == "no_coax")
    o = ModelCfg::Ctd::NO_COAX;
  else if (s == "all")
    o = ModelCfg::Ctd::ALL;
  else
    error("Invalid CTD option %s", s.c_str());
  return str;
}

}  // namespace mrna
