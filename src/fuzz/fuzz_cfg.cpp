// Copyright 2021 Eliot Courtney.
#include "fuzz/fuzz_cfg.h"

#include <fmt/core.h>

#include "compute/energy/energy_cfg.h"

namespace mrna::fuzz {

void RegisterOpts(ArgParse* args) {
  erg::RegisterOpts(args);
  args->RegisterOpt(OPT_FUZZ_BRUTE_MAX);
  args->RegisterOpt(OPT_FUZZ_MFE);
  args->RegisterOpt(OPT_FUZZ_MFE_RNASTRUCTURE);
  args->RegisterOpt(OPT_FUZZ_MFE_TABLE);
  args->RegisterOpt(OPT_FUZZ_SUBOPT);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_RNASTRUCTURE);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_STRUCS);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_DELTA);
  args->RegisterOpt(OPT_FUZZ_PARTITION);
  args->RegisterOpt(OPT_FUZZ_PARTITION_RNASTRUCTURE);
}

std::string FuzzCfg::Desc() const {
  std::string desc;
  desc += fmt::format("brute_max: {}\n", brute_max);
  desc += fmt::format("mfe: {}\n", mfe);
  desc += fmt::format("mfe_rnastructure: {}\n", mfe_rnastructure);
  desc += fmt::format("mfe_table: {}\n", mfe_table);
  desc += fmt::format("subopt: {}\n", subopt);
  desc += fmt::format("subopt_rnastructure: {}\n", subopt_rnastructure);
  desc += fmt::format("subopt_max: {}\n", subopt_strucs);
  desc += fmt::format("subopt_delta: {}\n", subopt_delta);
  desc += fmt::format("part: {}\n", part);
  desc += fmt::format("part_rnastructure: {}\n", part_rnastructure);
  return desc;
}

FuzzCfg FuzzCfg::FromArgParse(const ArgParse& args) {
  FuzzCfg cfg;
  args.MaybeSet(OPT_FUZZ_BRUTE_MAX, &cfg.brute_max);

  args.MaybeSet(OPT_FUZZ_MFE, &cfg.mfe);
  args.MaybeSet(OPT_FUZZ_MFE_RNASTRUCTURE, &cfg.mfe_rnastructure);
  args.MaybeSet(OPT_FUZZ_MFE_TABLE, &cfg.mfe_table);

  args.MaybeSet(OPT_FUZZ_SUBOPT, &cfg.subopt);
  args.MaybeSet(OPT_FUZZ_SUBOPT_RNASTRUCTURE, &cfg.subopt_rnastructure);
  args.MaybeSet(OPT_FUZZ_SUBOPT_STRUCS, &cfg.subopt_strucs);
  args.MaybeSet(OPT_FUZZ_SUBOPT_DELTA, &cfg.subopt_delta);

  args.MaybeSet(OPT_FUZZ_PARTITION, &cfg.part);
  args.MaybeSet(OPT_FUZZ_PARTITION_RNASTRUCTURE, &cfg.part_rnastructure);

  cfg.mfe = cfg.mfe || cfg.mfe_rnastructure;
  cfg.subopt = cfg.subopt || cfg.subopt_rnastructure;
  cfg.part = cfg.part || cfg.part_rnastructure;

  return cfg;
}

}  // namespace mrna::fuzz
