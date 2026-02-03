// Copyright 2022 Eliot Courtney.
#include "api/brute/brute_cfg.h"

#include <ostream>

#include "api/options.h"

namespace mrna::brute {

void RegisterOpts(ArgParse* args) {
  subopt::RegisterOpts(args);
  args->RegisterOpt(OPT_FOLD);
  args->RegisterOpt(OPT_SUBOPT);
  args->RegisterOpt(OPT_PFN);
}

BruteCfg BruteCfg::FromArgParse(const ArgParse& args) {
  BruteCfg cfg;
  args.MaybeSet(OPT_FOLD, &cfg.mfe);
  args.MaybeSet(OPT_SUBOPT, &cfg.subopt);
  args.MaybeSet(OPT_PFN, &cfg.pfn);
  cfg.subopt_cfg = subopt::SuboptCfg::FromArgParse(args);
  return cfg;
}

std::ostream& operator<<(std::ostream& str, const BruteCfg& o) {
  return str << "BruteCfg{mfe=" << o.mfe << ", subopt=" << o.subopt << ", pfn=" << o.pfn
             << ", subopt_cfg=" << o.subopt_cfg << "}";
}

}  // namespace mrna::brute
