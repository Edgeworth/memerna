// Copyright 2022 Eliot Courtney.
#include "compute/energy/energy_cfg.h"

#include "util/error.h"

namespace mrna::erg {

void RegisterOpts(ArgParse* args) {
  args->RegisterOpt(OPT_SEED);
  args->RegisterOpt(OPT_MEMERNA_DATA);
  args->RegisterOpt(OPT_ENERGY_MODEL);
  args->RegisterOpt(OPT_LONELY_PAIRS);
  args->RegisterOpt(OPT_BULGE_STATES);
  args->RegisterOpt(OPT_CTD);
}

std::string ModelPathFromArgParse(const ArgParse& args) {
  return ModelPathFromArgParse(args, args.Get(OPT_ENERGY_MODEL));
}

std::string ModelPathFromArgParse(const ArgParse& args, const std::string& model) {
  return ModelPath(args.Get(OPT_MEMERNA_DATA), model);
}

std::string ModelPath(const std::string& data_dir, const std::string& model) {
  return data_dir + "/model/" + model;
}

EnergyCfg EnergyCfg::FromArgParse(const ArgParse& args) {
  EnergyCfg cfg;
  args.MaybeSet(OPT_LONELY_PAIRS, &cfg.lonely_pairs);
  args.MaybeSet(OPT_CTD, &cfg.ctd);
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

}  // namespace mrna::erg
