// Copyright 2021 Eliot Courtney.
#include "fuzz/fuzz_cfg.h"

#include <fmt/core.h>

#include <string>

#include "api/bridge/bridge.h"
#include "api/ctx/backend_cfg.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "util/error.h"

namespace mrna::fuzz {

void RegisterOpts(ArgParse* args) {
  erg::RegisterOptsEnergyCfg(args);
  erg::RegisterOptsPseudofree(args);
  args->RegisterOpt(OPT_ENERGY_MODEL);
  args->RegisterOpt(OPT_ENERGY_PRECISION);
  args->RegisterOpt(OPT_MEMERNA_DATA);
  args->RegisterOpt(OPT_SEED);
  args->RegisterOpt(OPT_RAND_MIN_ENERGY);
  args->RegisterOpt(OPT_RAND_MAX_ENERGY);
  args->RegisterOpt(OPT_FUZZ_BRUTE_MAX);
  args->RegisterOpt(OPT_FUZZ_MFE);
  args->RegisterOpt(OPT_FUZZ_MFE_RNASTRUCTURE);
  args->RegisterOpt(OPT_FUZZ_MFE_TABLE);
  args->RegisterOpt(OPT_FUZZ_SUBOPT);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_RNASTRUCTURE);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_STRUCS);
  args->RegisterOpt(OPT_FUZZ_SUBOPT_DELTA);
  args->RegisterOpt(OPT_FUZZ_PFN);
  args->RegisterOpt(OPT_FUZZ_PFN_RNASTRUCTURE);
  args->RegisterOpt(OPT_FUZZ_PFN_SUBOPT);
  args->RegisterOpt(OPT_FUZZ_PFN_PQ_REL_EP);
  args->RegisterOpt(OPT_FUZZ_PFN_PQ_ABS_EP);
  args->RegisterOpt(OPT_FUZZ_PFN_PROB_REL_EP);
  args->RegisterOpt(OPT_FUZZ_PFN_PROB_ABS_EP);
  args->RegisterOpt(OPT_FUZZ_BACKENDS);
  args->RegisterOpt(OPT_FUZZ_RANDOM_SEEDS);
  args->RegisterOpt(OPT_FUZZ_RANDOM_PSEUDOFREE);
  args->RegisterOpt(OPT_FUZZ_RANDOM_PSEUDOFREE_MIN_ENERGY);
  args->RegisterOpt(OPT_FUZZ_RANDOM_PSEUDOFREE_MAX_ENERGY);
  args->RegisterOpt(mrna::bridge::OPT_RNASTRUCTURE_DATA);
}

RandomPseudofreeCfg::RandomPseudofreeCfg(Energy min_energy, Energy max_energy)
    : min_energy(min_energy), max_energy(max_energy) {
  verify(this->min_energy <= this->max_energy,
      "random pseudofree min energy {} > max energy {}", this->min_energy, this->max_energy);
}

RandomPseudofreeCfg RandomPseudofreeCfg::FromArgParse(const ArgParse& args) {
  const bool explicit_range = args.HasExplicit(OPT_FUZZ_RANDOM_PSEUDOFREE_MIN_ENERGY) ||
      args.HasExplicit(OPT_FUZZ_RANDOM_PSEUDOFREE_MAX_ENERGY);
  verify(!explicit_range || args.GetOr(OPT_FUZZ_RANDOM_PSEUDOFREE),
      "cannot set random pseudofree energy range without --random-pf");
  return {args.Get<Energy>(OPT_FUZZ_RANDOM_PSEUDOFREE_MIN_ENERGY),
      args.Get<Energy>(OPT_FUZZ_RANDOM_PSEUDOFREE_MAX_ENERGY)};
}

std::ostream& operator<<(std::ostream& str, const RandomPseudofreeCfg& o) {
  return str << "RandomPseudofreeCfg{"
             << "min_energy=" << o.min_energy << ", max_energy=" << o.max_energy << "}";
}

std::string FuzzCfg::Desc() const {
  std::string desc;
  desc += fmt::format("brute_max: {}\n", brute_max);
  desc += fmt::format("mfe: {}\n", mfe);
  desc += fmt::format("mfe_rnastructure: {}\n", mfe_rnastructure);
  desc += fmt::format("mfe_table: {}\n", mfe_table);
  desc += fmt::format("subopt: {}\n", subopt);
  desc += fmt::format("subopt_rnastructure: {}\n", subopt_rnastructure);
  desc += fmt::format("subopt_strucs: {}\n", subopt_strucs);
  desc += fmt::format("subopt_delta: {}\n", subopt_delta);
  desc += fmt::format("pfn: {}\n", pfn);
  desc += fmt::format("pfn_rnastructure: {}\n", pfn_rnastructure);
  desc += fmt::format("pfn_subopt: {}\n", pfn_subopt);
  desc += fmt::format("pfn_pq_rel_ep: {}\n", pfn_pq_rel_ep);
  desc += fmt::format("pfn_pq_abs_ep: {}\n", pfn_pq_abs_ep);
  desc += fmt::format("pfn_prob_rel_ep: {}\n", pfn_prob_rel_ep);
  desc += fmt::format("pfn_prob_abs_ep: {}\n", pfn_prob_abs_ep);
  desc += fmt::format("random_seeds: {}\n", random_seeds);
  desc += fmt::format("random_model_cfg: {}\n", random_model_cfg);
  desc += fmt::format("random_pseudofree: {}\n", random_pseudofree);
  desc += fmt::format("random_pseudofree_cfg: {}\n", random_pseudofree_cfg);
  desc += fmt::format("energy_cfg: {}\n", energy_cfg);
  desc += fmt::format("energy_model: {}\n", energy_model);
  desc += fmt::format("backend: ");
  for (const auto& backend : backends) desc += fmt::format("{},", backend);
  desc += "\n";
  desc += fmt::format("data_dir: {}\n", data_dir);
  desc += fmt::format("rnastructure_data_dir: {}\n", rnastructure_data_dir);
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

  args.MaybeSet(OPT_FUZZ_PFN, &cfg.pfn);
  args.MaybeSet(OPT_FUZZ_PFN_RNASTRUCTURE, &cfg.pfn_rnastructure);
  args.MaybeSet(OPT_FUZZ_PFN_SUBOPT, &cfg.pfn_subopt);
  args.MaybeSet(OPT_FUZZ_PFN_PQ_REL_EP, &cfg.pfn_pq_rel_ep);
  args.MaybeSet(OPT_FUZZ_PFN_PQ_ABS_EP, &cfg.pfn_pq_abs_ep);
  args.MaybeSet(OPT_FUZZ_PFN_PROB_REL_EP, &cfg.pfn_prob_rel_ep);
  args.MaybeSet(OPT_FUZZ_PFN_PROB_ABS_EP, &cfg.pfn_prob_abs_ep);

  cfg.mfe = cfg.mfe || cfg.mfe_rnastructure;
  cfg.subopt = cfg.subopt || cfg.subopt_rnastructure;
  cfg.pfn = cfg.pfn || cfg.pfn_rnastructure || cfg.pfn_subopt;

  args.MaybeSet(OPT_FUZZ_RANDOM_SEEDS, &cfg.random_seeds);
  args.MaybeSet(OPT_FUZZ_RANDOM_PSEUDOFREE, &cfg.random_pseudofree);
  cfg.random_model_cfg = RandomModelCfg::FromArgParse(args);
  verify(!(cfg.random_seeds && cfg.random_model_cfg.seed.has_value()),
      "cannot set fixed seed with random seeds");
  cfg.random_pseudofree_cfg = RandomPseudofreeCfg::FromArgParse(args);

  cfg.energy_cfg = erg::EnergyCfg::FromArgParse(args);
  cfg.energy_model = args.Get<erg::EnergyModelKind>(OPT_ENERGY_MODEL);
  cfg.backends = args.GetMultipleOr<BackendKind>(OPT_FUZZ_BACKENDS);
  cfg.pf_paired = args.GetMultipleOr<Energy>(erg::OPT_PAIRED_PSEUDOFREE);
  cfg.pf_unpaired = args.GetMultipleOr<Energy>(erg::OPT_UNPAIRED_PSEUDOFREE);
  verify(!cfg.random_pseudofree || (cfg.pf_paired.empty() && cfg.pf_unpaired.empty()),
      "cannot set pseudofree energies with random pseudofree");

  cfg.data_dir = args.Get<std::string>(OPT_MEMERNA_DATA);

#ifdef MRNA_USE_RNASTRUCTURE
  cfg.rnastructure_data_dir = args.Get<std::string>(bridge::OPT_RNASTRUCTURE_DATA);
#endif  // MRNA_USE_RNASTRUCTURE

  return cfg;
}

}  // namespace mrna::fuzz
