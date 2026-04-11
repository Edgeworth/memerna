// Copyright 2024 Eliot Courtney.
#include "api/ctx/backend_cfg.h"

#include <fmt/format.h>

#include <ostream>
#include <string>

#include "util/util.h"

namespace mrna {

RandomModelCfg::RandomModelCfg(
    std::optional<uint_fast32_t> seed, Energy min_energy, Energy max_energy)
    : seed(seed), min_energy(min_energy), max_energy(max_energy) {
  verify(this->min_energy <= this->max_energy, "random model min energy {} > max energy {}",
      this->min_energy, this->max_energy);
}

RandomModelCfg RandomModelCfg::FromArgParse(const ArgParse& args) {
  return {args.MaybeGet<uint_fast32_t>(OPT_SEED),
      args.Get<Energy>(OPT_RAND_MIN_ENERGY),
      args.Get<Energy>(OPT_RAND_MAX_ENERGY)};
}

std::ostream& operator<<(std::ostream& str, const RandomModelCfg& o) {
  str << "RandomModelCfg{";
  if (o.seed.has_value()) str << "seed=" << *o.seed << ", ";
  return str << "min_energy=" << o.min_energy << ", max_energy=" << o.max_energy << "}";
}

BackendCfg BackendCfg::FromArgParse(const ArgParse& args) {
  BackendCfg cfg{
      .energy_model = args.Get<erg::EnergyModelKind>(OPT_ENERGY_MODEL),
      .precision = args.Get<int>(OPT_ENERGY_PRECISION),
  };
  const auto random_cfg = RandomModelCfg::FromArgParse(args);
  if (random_cfg.seed.has_value()) {
    cfg.data_src = random_cfg;
  } else {
    verify(!args.HasExplicit(OPT_RAND_MIN_ENERGY) && !args.HasExplicit(OPT_RAND_MAX_ENERGY),
        "cannot set random model energy range without --seed");
    cfg.data_src = args.Get<std::string>(OPT_MEMERNA_DATA);
  }
  verify(cfg.precision == MRNA_ENERGY_PRECISION, "unsupported energy precision: {}, built with {}",
      cfg.precision, MRNA_ENERGY_PRECISION);
  return cfg;
}

std::optional<std::string> BackendCfg::ModelPath(BackendKind backend) const {
  return std::visit(
      overloaded{
          [](std::monostate) -> std::optional<std::string> { return std::nullopt; },
          [this, backend](const std::string& data_dir) -> std::optional<std::string> {
            return fmt::format("{}/model/{}-p{}-{}", data_dir, energy_model, precision, backend);
          },
          [](const RandomModelCfg&) -> std::optional<std::string> { return std::nullopt; },
      },
      data_src);
}

void RegisterOptsBackendCfg(ArgParse* args) {
  erg::RegisterOptsEnergyCfg(args);
  args->RegisterOpt(OPT_ENERGY_MODEL);
  args->RegisterOpt(OPT_ENERGY_PRECISION);
  args->RegisterOpt(OPT_BACKEND);
  args->RegisterOpt(OPT_MEMERNA_DATA);
  args->RegisterOpt(OPT_SEED);
  args->RegisterOpt(OPT_RAND_MIN_ENERGY);
  args->RegisterOpt(OPT_RAND_MAX_ENERGY);
}

}  // namespace mrna
