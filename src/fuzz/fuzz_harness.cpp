// Copyright 2022 Eliot Courtney.
#include "fuzz/fuzz_harness.h"

#include <fmt/core.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "api/ctx/algorithm.h"
#include "api/ctx/backend.h"
#include "api/energy/pseudofree_cfg.h"
#include "fuzz/fuzz_cfg.h"
#include "model/primary.h"
#include "util/log.h"

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

FuzzHarness::FuzzHarness(FuzzCfg fuzz_cfg, bool should_log)
    : fuzz_cfg_(std::move(fuzz_cfg)), e_(std::random_device{}()), should_log_(should_log) {
#ifdef MRNA_USE_RNASTRUCTURE
  rstr_ =
      std::make_shared<bridge::RNAstructure>(fuzz_cfg_.rnastructure_data_dir, /*use_lyngso=*/false);
#endif  // MRNA_USE_RNASTRUCTURE
  if (should_log_) fmt::print("Fuzzing with config: {}\n", fuzz_cfg_.Desc());

  backend_cfg_ = BackendCfg{
      .energy_model = fuzz_cfg_.energy_model,
      .precision = MRNA_ENERGY_PRECISION,
      .data_src = fuzz_cfg_.data_dir,
  };
}

Error FuzzHarness::Run(const Primary& r, const erg::PseudofreeCfg& pf) {
  MaybeLoadBackends(pf);
  auto groups = EquivalenceClasses(ms_);
  if (should_log_) {
    loginfo("backend equivalence classes: {}", groups.size());
    for (int i = 0; i < static_cast<int>(groups.size()); ++i)
      loginfo("equivalence class {}: {}", i, DescribeEquivalenceClass(groups[i]));
  }

#ifdef MRNA_USE_RNASTRUCTURE
  const bool rnastructure =
      fuzz_cfg_.mfe_rnastructure || fuzz_cfg_.subopt_rnastructure || fuzz_cfg_.pfn_rnastructure;
  verify(!rnastructure || groups.size() == 1,
      "RNAstructure comparison requires exactly 1 backend equivalence class, got {}",
      groups.size());
#endif  // MRNA_USE_RNASTRUCTURE

  Error errors;
  for (const auto& group : groups) {
    FuzzInvocation invoc(r, group, backend_cfg_, pf, fuzz_cfg_, should_log_);
#ifdef MRNA_USE_RNASTRUCTURE
    invoc.set_rnastructure(rstr_);
#endif  // MRNA_USE_RNASTRUCTURE
    auto local = invoc.Run();
    errors.insert(
        errors.end(), std::make_move_iterator(local.begin()), std::make_move_iterator(local.end()));
  }

  should_log_ = false;
  return errors;
}

void FuzzHarness::MaybeLoadBackends(const erg::PseudofreeCfg& pf) {
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

  if (fuzz_cfg_.backends.empty()) {
    for (const auto& backend : EnumValues<BackendKind>()) {
      std::string reason;
      if (!BackendIsSupported(backend, backend_cfg_, fuzz_cfg_.energy_cfg, pf, &reason)) {
        if (should_log_) logwarn("backend NOT loaded: {}: {}", backend, reason);
        continue;
      }
      ms_.push_back(BackendFromBackendCfg(backend, backend_cfg_));
    }
  } else {
    for (const auto& backend : fuzz_cfg_.backends) {
      ms_.push_back(BackendFromBackendCfg(backend, backend_cfg_));
    }
  }
}

}  // namespace mrna::fuzz
