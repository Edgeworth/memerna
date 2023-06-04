// Copyright 2022 Eliot Courtney.
#include "fuzz/fuzz_harness.h"

#include <fmt/core.h>

#include <cinttypes>
#include <ctime>
#include <memory>
#include <string>

#include "api/energy/model.h"
#include "fuzz/fuzz_cfg.h"
#include "model/primary.h"
#include "util/string.h"

namespace mrna::fuzz {

FuzzHarness::FuzzHarness(FuzzCfg cfg) : cfg_(cfg), e_(uint_fast32_t(time(nullptr))) {
#ifdef USE_RNASTRUCTURE
  rstr_ = std::make_shared<bridge::RNAstructure>(cfg.rnastructure_data_dir, false);
#endif  // USE_RNASTRUCTURE
  fmt::print("Fuzzing with config: {}\n", cfg_.Desc());
}

FuzzInvocation FuzzHarness::CreateInvocation(const Primary& r) {
  uint_fast32_t seed = e_();
  MaybeLoadModels(seed);

  FuzzInvocation invoc(r, ems_, cfg_, 0);
#ifdef USE_RNASTRUCTURE
  invoc.set_rnastructure(rstr_);
#endif  // USE_RNASTRUCTURE
  return invoc;
}

void FuzzHarness::MaybeLoadModels(uint_fast32_t seed) {
  // Don't reload if already loaded and not randomising.
  if (!ems_.empty() && !cfg_.random_models) return;
  ems_.clear();

  for (const auto& model_name : cfg_.model_names) {
    ems_.push_back(LoadModel(model_name, seed));
  }
}

erg::EnergyModelPtr FuzzHarness::LoadModel(
    const std::string& model_name, uint_fast32_t seed) const {
  auto kind = Conv<erg::ModelKind>(model_name);
  if (cfg_.seed >= 0) return erg::Random(kind, cfg_.seed);
  if (cfg_.random_models) return erg::Random(kind, seed);
  return erg::FromDir(cfg_.data_dir, model_name);
}

}  // namespace mrna::fuzz
