// Copyright 2022 Eliot Courtney.
#include "fuzz/fuzz_harness.h"

#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <memory>
#include <utility>
#include <vector>

#include "api/ctx/backend.h"
#include "fuzz/fuzz_cfg.h"
#include "model/primary.h"

namespace mrna::fuzz {

namespace {

std::vector<std::vector<BackendModelPtr>> EquivalenceClasses(
    const std::vector<BackendModelPtr>& ms) {
  std::vector<std::vector<BackendModelPtr>> groups;
  for (const auto& m : ms) {
    bool found = false;
    for (auto& group : groups) {
      if (IsEquivalent(group[0], m)) {
        group.push_back(m);
        found = true;
        break;
      }
    }
    if (!found) groups.push_back({m});
  }
  return groups;
}

std::string DescribeEquivalenceClass(const std::vector<BackendModelPtr>& group) {
  std::string desc;
  for (int i = 0; i < static_cast<int>(group.size()); ++i) {
    if (i > 0) desc += ", ";
    desc += fmt::format("{}", GetBackendKind(group[i]));
  }
  return desc;
}

}  // namespace

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

Error FuzzHarness::Run(const Primary& r, erg::PseudofreeCfg pf) {
  MaybeLoadBackends();
  auto groups = EquivalenceClasses(ms_);
  if (first_invocation_) {
    spdlog::info("backend equivalence classes: {}", groups.size());
    for (int i = 0; i < static_cast<int>(groups.size()); ++i)
      spdlog::info("equivalence class {}: {}", i, DescribeEquivalenceClass(groups[i]));
  }

#ifdef USE_RNASTRUCTURE
  const bool rnastructure =
      fuzz_cfg_.mfe_rnastructure || fuzz_cfg_.subopt_rnastructure || fuzz_cfg_.pfn_rnastructure;
  verify(
      !rnastructure || groups.size() == 1,
      "RNAstructure comparison requires exactly 1 backend equivalence class, got {}",
      groups.size());
#endif  // USE_RNASTRUCTURE

  Error errors;
  for (const auto& group : groups) {
    FuzzInvocation invoc(r, group, backend_cfg_, pf, fuzz_cfg_, first_invocation_);
#ifdef USE_RNASTRUCTURE
    invoc.set_rnastructure(rstr_);
#endif  // USE_RNASTRUCTURE
    auto local = invoc.Run();
    errors.insert(errors.end(), std::make_move_iterator(local.begin()), std::make_move_iterator(local.end()));
  }

  first_invocation_ = false;
  return errors;
}

void FuzzHarness::MaybeLoadBackends() {
  // Don't reload if already loaded and not randomising.
  if (!ms_.empty() && !fuzz_cfg_.random_models) return;
  ms_.clear();

  if (fuzz_cfg_.seed.has_value()) {
    backend_cfg_.data_src = *fuzz_cfg_.seed;
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
