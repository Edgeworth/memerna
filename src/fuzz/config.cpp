// Copyright 2021 E.
#include "fuzz/config.h"

#include <cstdlib>

#include "util/error.h"
#include "util/string.h"

namespace mrna::fuzz {

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
  cfg.random_model = args.HasFlag("random");
  cfg.table_check = args.HasFlag("table-check");
  if (args.HasFlag("brute-cutoff")) cfg.brute_cutoff = atoi(args.GetOption("brute-cutoff").c_str());
  if (args.HasFlag("brute-subopt-max"))
    cfg.brute_subopt_max = atoi(args.GetOption("brute-subopt-max").c_str());

  cfg.mfe_rnastructure = args.HasFlag("mfe-rnastructure");

  cfg.subopt = args.HasFlag("subopt");
  cfg.subopt_rnastructure = args.HasFlag("subopt-rnastructure");
  if (args.HasFlag("subopt-max")) cfg.subopt_max = atoi(args.GetOption("subopt-max").c_str());
  if (args.HasFlag("subopt-delta")) cfg.subopt_delta = atoi(args.GetOption("subopt-delta").c_str());

  cfg.partition = args.HasFlag("partition");
  cfg.partition_rnastructure = args.HasFlag("partition-rnastructure");

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
