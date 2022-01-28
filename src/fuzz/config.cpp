// Copyright 2021 Eliot Courtney.
#include "fuzz/config.h"

#include <cstdlib>

#include "compute/energy/energy.h"
#include "util/error.h"
#include "util/string.h"

namespace mrna::fuzz {

void RegisterOpts(ArgParse* args) {
  energy::RegisterOpts(args);
  args->RegisterOpt(OPT_RANDOM);
  args->RegisterOpt(OPT_TABLE_CHECK);
  args->RegisterOpt(OPT_BRUTE_CUTOFF);
  args->RegisterOpt(OPT_BRUTE_SUBOPT_MAX);
  args->RegisterOpt(OPT_MFE_RNASTRUCTURE);
  args->RegisterOpt(OPT_SUBOPT);
  args->RegisterOpt(OPT_SUBOPT_RNASTRUCTURE);
  args->RegisterOpt(OPT_SUBOPT_MAX);
  args->RegisterOpt(OPT_SUBOPT_DELTA);
  args->RegisterOpt(OPT_PARTITION);
  args->RegisterOpt(OPT_PARTITION_RNASTRUCTURE);
}

std::string FuzzCfg::Describe() {
  std::string desc;
  if (random_model)
    desc += "random energy models";
  else
    desc += "specified energy model";
  if (subopt) {
    desc += " - testing suboptimal";
    if (subopt_rnastructure) desc += " (including rnastructure)";
    desc += sfmt(" max: %d delta: %d", subopt_max, subopt_delta);
  }
  desc += sfmt(" - brute cutoff: %d subopt-max: %d", brute_cutoff, brute_subopt_max);
  if (partition) {
    desc += " - testing partition";
    if (partition_rnastructure) desc += " (including rnastructure)";
  }
  return desc;
}

FuzzCfg FuzzCfg::FromArgParse(const ArgParse& args) {
  FuzzCfg cfg;
  cfg.random_model = args.Has(OPT_RANDOM);
  cfg.table_check = args.Has(OPT_TABLE_CHECK);
  if (args.Has(OPT_BRUTE_CUTOFF)) cfg.brute_cutoff = atoi(args.Get(OPT_BRUTE_CUTOFF).c_str());
  if (args.Has(OPT_BRUTE_SUBOPT_MAX))
    cfg.brute_subopt_max = atoi(args.Get(OPT_BRUTE_SUBOPT_MAX).c_str());

  cfg.mfe_rnastructure = args.Has(OPT_MFE_RNASTRUCTURE);

  cfg.subopt = args.Has(OPT_SUBOPT);
  cfg.subopt_rnastructure = args.Has(OPT_SUBOPT_RNASTRUCTURE);
  if (args.Has(OPT_SUBOPT_MAX)) cfg.subopt_max = atoi(args.Get(OPT_SUBOPT_MAX).c_str());
  if (args.Has(OPT_SUBOPT_DELTA)) cfg.subopt_delta = atoi(args.Get(OPT_SUBOPT_DELTA).c_str());

  cfg.partition = args.Has(OPT_PARTITION);
  cfg.partition_rnastructure = args.Has(OPT_PARTITION_RNASTRUCTURE);

  verify(!cfg.subopt_rnastructure || cfg.subopt,
      "suboptimal folding testing must be enabled to test rnastructure suboptimal folding");
  verify(!(cfg.random_model && cfg.subopt_rnastructure),
      "cannot use a random energy model with rnastructure");
  verify(cfg.mfe_rnastructure || (!cfg.subopt_rnastructure && !cfg.partition_rnastructure),
      "rnastructure must be enabled to use it for suboptimal or partition");
  verify(!cfg.random_model || !cfg.mfe_rnastructure,
      "rnastructure testing does not support random models");
  return cfg;
}

}  // namespace mrna::fuzz
