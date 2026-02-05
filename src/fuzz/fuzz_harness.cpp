// Copyright 2022 Eliot Courtney.
#include "fuzz/fuzz_harness.h"

#include <fmt/core.h>

#include <memory>
#include <utility>
#include <vector>

#include "api/ctx/backend.h"
#include "fuzz/fuzz_cfg.h"
#include "model/primary.h"

namespace mrna::fuzz {

FuzzHarness::FuzzHarness(FuzzCfg fuzz_cfg)
    : fuzz_cfg_(std::move(fuzz_cfg)), e_(std::random_device{}()) {
#ifdef USE_RNASTRUCTURE
  rstr_ =
      std::make_shared<bridge::RNAstructure>(fuzz_cfg_.rnastructure_data_dir, /*use_lyngso=*/false);
#endif  // USE_RNASTRUCTURE
  fmt::print("Fuzzing with config: {}\n", fuzz_cfg_.Desc());

  backend_cfg_ = BackendCfg{
      .energy_model = fuzz_cfg_.energy_model,
      .precision = ENERGY_PRECISION,
      .data_src = fuzz_cfg_.data_dir,
  };
}

FuzzInvocation FuzzHarness::CreateInvocation(const Primary& r, erg::PseudofreeCfg pf) {
  MaybeLoadBackends();

  FuzzInvocation invoc(r, ms_, backend_cfg_, std::move(pf), fuzz_cfg_);
#ifdef USE_RNASTRUCTURE
  invoc.set_rnastructure(rstr_);
#endif  // USE_RNASTRUCTURE
  return invoc;
}

void FuzzHarness::MaybeLoadBackends() {
  // Don't reload if already loaded and not randomising.
  if (!ms_.empty() && !fuzz_cfg_.random_models) return;
  ms_.clear();

  if (fuzz_cfg_.seed >= 0) {
    backend_cfg_.data_src = static_cast<uint_fast32_t>(fuzz_cfg_.seed);
  } else if (fuzz_cfg_.random_models) {
    backend_cfg_.data_src = static_cast<uint_fast32_t>(e_());
  } else {
    backend_cfg_.data_src = fuzz_cfg_.data_dir;
  }

  for (const auto& backend : fuzz_cfg_.backends) {
    ms_.push_back(BackendFromBackendCfg(backend, backend_cfg_));
  }
}

}  // namespace mrna::fuzz
