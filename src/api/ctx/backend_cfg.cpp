// Copyright 2024 Eliot Courtney.
#include "api/ctx/backend_cfg.h"

#include <fmt/format.h>

namespace mrna {

BackendCfg BackendCfg::FromArgParse(const ArgParse& args) {
  BackendCfg cfg{
      .energy_model = args.Get<erg::EnergyModelKind>(OPT_ENERGY_MODEL),
      .precision = args.Get<int>(OPT_ENERGY_PRECISION),
      .backend = args.Get<BackendKind>(OPT_BACKEND),
      .data_dir = args.Get<std::string>(OPT_MEMERNA_DATA),
      .seed = args.MaybeGet<uint_fast32_t>(OPT_SEED),
      .energy_cfg = erg::EnergyCfg::FromArgParse(args),
      .pf_paired = args.GetMultipleOr<Energy>(OPT_PAIRED_PSEUDOFREE),
      .pf_unpaired = args.GetMultipleOr<Energy>(OPT_UNPAIRED_PSEUDOFREE),
  };
  verify(cfg.precision == ENERGY_PRECISION, "unsupported energy precision: {}, built with {}",
      cfg.precision, ENERGY_PRECISION);
  return cfg;
}

std::string BackendCfg::BackendDataPath() const {
  return fmt::format("{}/model/{}-p{}-{}", data_dir, energy_model, precision, backend);
}

void RegisterOptsBackendCfg(ArgParse* args) {
  erg::RegisterOptsEnergyCfg(args);
  args->RegisterOpt(OPT_ENERGY_MODEL);
  args->RegisterOpt(OPT_ENERGY_PRECISION);
  args->RegisterOpt(OPT_BACKEND);
  args->RegisterOpt(OPT_MEMERNA_DATA);
  args->RegisterOpt(OPT_SEED);
  args->RegisterOpt(OPT_PAIRED_PSEUDOFREE);
  args->RegisterOpt(OPT_UNPAIRED_PSEUDOFREE);
}

}  // namespace mrna
