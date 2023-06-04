// Copyright 2022 Eliot Courtney.
#include "api/energy/energy_cfg.h"

#include <ostream>

#include "util/error.h"

namespace mrna::erg {

void RegisterOptsEnergyCfg(ArgParse* args) {
  args->RegisterOpt(OPT_LONELY_PAIRS);
  args->RegisterOpt(OPT_BULGE_STATES);
  args->RegisterOpt(OPT_CTD);
}

EnergyCfg EnergyCfg::FromArgParse(const ArgParse& args) {
  EnergyCfg cfg;
  args.MaybeSet(OPT_LONELY_PAIRS, &cfg.lonely_pairs);
  args.MaybeSet(OPT_CTD, &cfg.ctd);
  return cfg;
}

std::ostream& operator<<(std::ostream& str, const EnergyCfg& o) {
  return str << "EnergyCfg{"
             << "lonely_pairs=" << o.lonely_pairs << ", ctd=" << o.ctd
             << ", bulge_states=" << o.bulge_states << "}";
}

std::ostream& operator<<(std::ostream& str, const EnergyCfg::LonelyPairs& o) {
  switch (o) {
  case EnergyCfg::LonelyPairs::OFF: return str << "off";
  case EnergyCfg::LonelyPairs::HEURISTIC: return str << "heuristic";
  case EnergyCfg::LonelyPairs::ON: return str << "on";
  default: bug();
  }
}

std::ostream& operator<<(std::ostream& str, const EnergyCfg::Ctd& o) {
  switch (o) {
  case EnergyCfg::Ctd::NONE: return str << "none";
  case EnergyCfg::Ctd::NO_COAX: return str << "no-coax";
  case EnergyCfg::Ctd::ALL: return str << "all";
  default: bug();
  }
}

std::istream& operator>>(std::istream& str, EnergyCfg::LonelyPairs& o) {
  std::string s;
  str >> s;
  if (s == "off")
    o = EnergyCfg::LonelyPairs::OFF;
  else if (s == "heuristic")
    o = EnergyCfg::LonelyPairs::HEURISTIC;
  else if (s == "on")
    o = EnergyCfg::LonelyPairs::ON;
  else
    fatal("Invalid CTD option {}", s);
  return str;
}

std::istream& operator>>(std::istream& str, EnergyCfg::Ctd& o) {
  std::string s;
  str >> s;
  if (s == "none")
    o = EnergyCfg::Ctd::NONE;
  else if (s == "no-coax")
    o = EnergyCfg::Ctd::NO_COAX;
  else if (s == "all")
    o = EnergyCfg::Ctd::ALL;
  else
    fatal("Invalid CTD option {}", s);
  return str;
}

}  // namespace mrna::erg
