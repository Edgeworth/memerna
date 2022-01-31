// Copyright 2021 E.
#include "fuzz/config.h"

#include "compute/energy/config.h"
#include "util/string.h"

namespace mrna::fuzz {

void RegisterOpts(ArgParse* args) {
  energy::RegisterOpts(args);
  args->RegisterOpt(OPT_FUZZ_BRUTE_MAX);
  args->RegisterOpt(OPT_FUZZ_MFE);
  args->RegisterOpt(OPT_FUZZ_MFE_RNASTRUCTURE);
  args->RegisterOpt(OPT_FUZZ_MFE_BRUTE);
  args->RegisterOpt(OPT_FUZZ_MFE_TABLE);
  args->RegisterOpt(OPT_FUZZ_SUBOPT);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_RNASTRUCTURE);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_BRUTE);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_MAX);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_DELTA);
  args->RegisterOpt(OPT_FUZZ_PARTITION);
  args->RegisterOpt(OPT_FUZZ_PARTITION_RNASTRUCTURE);
  args->RegisterOpt(OPT_FUZZ_PARTITION_BRUTE);
}

std::string FuzzCfg::Desc() {
  std::string desc;
  desc += sfmt("brute_max: %d\n", brute_max);
  desc += sfmt("mfe: %d\n", mfe);
  desc += sfmt("mfe_rnastructure: %d\n", mfe_rnastructure);
  desc += sfmt("mfe_brute: %d\n", mfe_brute);
  desc += sfmt("mfe_table: %d\n", mfe_table);
  desc += sfmt("subopt: %d\n", subopt);
  desc += sfmt("subopt_rnastructure: %d\n", subopt_rnastructure);
  desc += sfmt("subopt_brute: %d\n", subopt_brute);
  desc += sfmt("subopt_max: %d\n", subopt_max);
  desc += sfmt("subopt_delta: %d\n", subopt_delta);
  desc += sfmt("part: %d\n", part);
  desc += sfmt("part_rnastructure: %d\n", part_rnastructure);
  desc += sfmt("part_brute: %d\n", part_brute);
  return desc;
}

FuzzCfg FuzzCfg::FromArgParse(const ArgParse& args) {
  FuzzCfg cfg;
  args.MaybeSet(OPT_FUZZ_BRUTE_MAX, &cfg.brute_max);

  args.MaybeSet(OPT_FUZZ_MFE, &cfg.mfe);
  args.MaybeSet(OPT_FUZZ_MFE_RNASTRUCTURE, &cfg.mfe_rnastructure);
  args.MaybeSet(OPT_FUZZ_MFE_BRUTE, &cfg.mfe_brute);
  args.MaybeSet(OPT_FUZZ_MFE_TABLE, &cfg.mfe_table);

  args.MaybeSet(OPT_FUZZ_SUBOPT, &cfg.subopt);
  args.MaybeSet(OPT_FUZZ_SUBOPT_RNASTRUCTURE, &cfg.subopt_rnastructure);
  args.MaybeSet(OPT_FUZZ_SUBOPT_BRUTE, &cfg.subopt_brute);
  args.MaybeSet(OPT_FUZZ_SUBOPT_MAX, &cfg.subopt_max);
  args.MaybeSet(OPT_FUZZ_SUBOPT_DELTA, &cfg.subopt_delta);

  args.MaybeSet(OPT_FUZZ_PARTITION, &cfg.part);
  args.MaybeSet(OPT_FUZZ_PARTITION_RNASTRUCTURE, &cfg.part_rnastructure);
  args.MaybeSet(OPT_FUZZ_PARTITION_BRUTE, &cfg.part_brute);

  cfg.mfe = cfg.mfe || cfg.mfe_rnastructure || cfg.mfe_brute;
  cfg.subopt = cfg.subopt || cfg.subopt_rnastructure || cfg.subopt_brute;
  cfg.part = cfg.part || cfg.part_rnastructure || cfg.part_brute;

  return cfg;
}

}  // namespace mrna::fuzz
