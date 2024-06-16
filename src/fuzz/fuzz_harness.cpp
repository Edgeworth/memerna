// Copyright 2022 Eliot Courtney.
#include "fuzz/fuzz_harness.h"

#include <fmt/core.h>

#include <memory>

#include "api/ctx/backend.h"
#include "fuzz/fuzz_cfg.h"
#include "model/primary.h"

namespace mrna::fuzz {

FuzzHarness::FuzzHarness(FuzzCfg fuzz_cfg)
    : fuzz_cfg_(std::move(fuzz_cfg)), e_(std::random_device{}()) {
#ifdef USE_RNASTRUCTURE
  rstr_ = std::make_shared<bridge::RNAstructure>(fuzz_cfg_.rnastructure_data_dir, false);
#endif  // USE_RNASTRUCTURE
  fmt::print("Fuzzing with config: {}\n", fuzz_cfg_.Desc());

  backend_cfg_ = BackendCfg{
      .energy_model = fuzz_cfg_.energy_model,
      .data_dir = fuzz_cfg_.data_dir,
      .seed = 0,
      .energy_cfg = fuzz_cfg_.energy_cfg,
  };
}

FuzzInvocation FuzzHarness::CreateInvocation(
    const Primary& r, std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired) {
  MaybeLoadBackends(std::move(pf_paired), std::move(pf_unpaired));

  FuzzInvocation invoc(r, ms_, fuzz_cfg_);
#ifdef USE_RNASTRUCTURE
  invoc.set_rnastructure(rstr_);
#endif  // USE_RNASTRUCTURE
  return invoc;
}

void FuzzHarness::MaybeLoadBackends(
    std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired) {
  // Don't reload if already loaded and not randomising.
  if (!ms_.empty() && !fuzz_cfg_.random_models) {
    bool pf_paired_changed = backend_cfg_.pf_paired != pf_paired;
    bool pf_unpaired_changed = backend_cfg_.pf_unpaired != pf_unpaired;
    if (!pf_paired_changed && !pf_unpaired_changed) return;
  }
  ms_.clear();

  backend_cfg_.seed = std::nullopt;
  if (fuzz_cfg_.seed >= 0) backend_cfg_.seed = fuzz_cfg_.seed;
  if (fuzz_cfg_.random_models) backend_cfg_.seed = e_();

  backend_cfg_.pf_paired = std::move(pf_paired);
  backend_cfg_.pf_unpaired = std::move(pf_unpaired);

  for (const auto& backend : fuzz_cfg_.backends) {
    backend_cfg_.backend = backend;
    ms_.push_back(BackendFromBackendCfg(backend_cfg_));
  }
}

}  // namespace mrna::fuzz
