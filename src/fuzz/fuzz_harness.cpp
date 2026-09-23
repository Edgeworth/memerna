// Copyright 2022 Eliot Courtney.
#include "fuzz/fuzz_harness.h"

#include <fmt/format.h>

#include <cstdint>
#include <iterator>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "api/ctx/algorithm.h"
#include "api/ctx/backend.h"
#include "api/ctx/backend_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_invocation.h"
#include "model/primary.h"
#include "util/enum.h"
#include "util/log.h"

#ifdef MRNA_USE_RNASTRUCTURE
#include <memory>

#include "api/bridge/rnastructure.h"
#include "util/error.h"
#endif  // MRNA_USE_RNASTRUCTURE

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
  fuzz_cfg_.Verify();
#ifdef MRNA_USE_RNASTRUCTURE
  if (fuzz_cfg_.UsesRNAstructure())
    rstr_ = std::make_shared<bridge::RNAstructure>(
        fuzz_cfg_.rnastructure_data_dir, /*use_lyngso=*/false);
#endif  // MRNA_USE_RNASTRUCTURE
  if (should_log_) loginfo("Fuzzing with config: {}", fuzz_cfg_.Desc());

  backend_cfg_ = BackendCfg{
      .energy_model = fuzz_cfg_.energy_model,
      .precision = MRNA_ENERGY_PRECISION,
      .data_src = fuzz_cfg_.data_dir,
  };
}

std::optional<Error> FuzzHarness::Run(const Primary& r, const erg::PseudofreeCfg& pf) {
  MaybeLoadBackends(pf);
  if (ms_.empty()) return std::nullopt;
  auto groups = EquivalenceClasses(ms_);
  if (should_log_) {
    loginfo("backend equivalence classes: {}", groups.size());
    for (int i = 0; i < static_cast<int>(groups.size()); ++i)
      loginfo("equivalence class {}: {}", i, DescribeEquivalenceClass(groups[i]));
  }

#ifdef MRNA_USE_RNASTRUCTURE
  verify(!fuzz_cfg_.UsesRNAstructure() || groups.size() == 1,
      "RNAstructure comparison requires exactly 1 backend equivalence class, got {}",
      groups.size());
#endif  // MRNA_USE_RNASTRUCTURE

  Error errors;
  for (const auto& group : groups) {
    FuzzInvocation invoc(r, group, backend_cfg_, pf, fuzz_cfg_, should_log_);
#ifdef MRNA_USE_RNASTRUCTURE
    if (fuzz_cfg_.UsesRNAstructure()) invoc.set_rnastructure(rstr_);
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
  if (!ms_.empty() && !fuzz_cfg_.random_seeds) return;
  ms_.clear();

  if (fuzz_cfg_.random_seeds) {
    backend_cfg_.data_src = RandomModelCfg(static_cast<uint_fast32_t>(e_()),
        fuzz_cfg_.random_model_cfg.min_energy, fuzz_cfg_.random_model_cfg.max_energy);
  } else if (fuzz_cfg_.random_model_cfg.seed.has_value()) {
    backend_cfg_.data_src = fuzz_cfg_.random_model_cfg;
  } else {
    backend_cfg_.data_src = fuzz_cfg_.data_dir;
  }

  const auto backends = fuzz_cfg_.backends.empty() ? EnumValues<BackendKind>() : fuzz_cfg_.backends;
  for (const auto& backend : backends) {
    std::string reason;
    if (!BackendIsSupported(backend, backend_cfg_, fuzz_cfg_.energy_cfg, pf, &reason)) {
      if (should_log_) logwarn("backend NOT loaded: {}: {}", backend, reason);
      continue;
    }
    ms_.push_back(BackendFromBackendCfg(backend, backend_cfg_));
  }
}

}  // namespace mrna::fuzz
