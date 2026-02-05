// Copyright 2024 Eliot Courtney.
#include "api/ctx/backend_cfg.h"

#include <fmt/format.h>

#include <string>

#include "util/util.h"

namespace mrna {

BackendCfg BackendCfg::FromArgParse(const ArgParse& args) {
  BackendCfg cfg{
      .energy_model = args.Get<erg::EnergyModelKind>(OPT_ENERGY_MODEL),
      .precision = args.Get<int>(OPT_ENERGY_PRECISION),
      .data_src = args.Has(OPT_SEED)
          ? std::variant<std::monostate, std::string, uint_fast32_t>{args.Get<uint_fast32_t>(
                OPT_SEED)}
          : std::variant<std::monostate, std::string, uint_fast32_t>{args.Get<std::string>(
                OPT_MEMERNA_DATA)},
  };
  verify(cfg.precision == ENERGY_PRECISION, "unsupported energy precision: {}, built with {}",
      cfg.precision, ENERGY_PRECISION);
  return cfg;
}

std::optional<std::string> BackendCfg::ModelPath(BackendKind backend) const {
  return std::visit(overloaded{
                        [](std::monostate) -> std::optional<std::string> { return std::nullopt; },
                        [this, backend](const std::string& data_dir) -> std::optional<std::string> {
                          return fmt::format(
                              "{}/model/{}-p{}-{}", data_dir, energy_model, precision, backend);
                        },
                        [](uint_fast32_t) -> std::optional<std::string> { return std::nullopt; },
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
}

}  // namespace mrna
