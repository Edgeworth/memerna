// Copyright 2022 Eliot Courtney.
#ifndef FUZZ_FUZZ_HARNESS_H_
#define FUZZ_FUZZ_HARNESS_H_

#include <cstdint>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <variant>
#include <vector>

#include "api/energy/pseudofree_cfg.h"
#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_invocation.h"
#include "model/primary.h"
#include "util/util.h"

#ifdef USE_RNASTRUCTURE
#include "api/bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::fuzz {

class FuzzHarness {
 public:
  explicit FuzzHarness(FuzzCfg fuzz_cfg);

  FuzzInvocation CreateInvocation(const Primary& r, erg::PseudofreeCfg pf);
  std::mt19937& e() { return e_; }

  [[nodiscard]] constexpr std::optional<uint_fast32_t> last_seed() const {
    return std::visit(
        overloaded{
            [](std::monostate) -> std::optional<uint_fast32_t> { return std::nullopt; },
            [](const std::string&) -> std::optional<uint_fast32_t> { return std::nullopt; },
            [](uint_fast32_t seed) -> std::optional<uint_fast32_t> { return seed; },
        },
        backend_cfg_.data_src);
  }

 private:
  std::vector<BackendModelPtr> ms_;
  BackendCfg backend_cfg_;
  FuzzCfg fuzz_cfg_;
  std::mt19937 e_;
  bool first_invocation_ = true;

#ifdef USE_RNASTRUCTURE
  std::shared_ptr<bridge::RNAstructure> rstr_;
#endif  // USE_RNASTRUCTURE

  void MaybeLoadBackends();
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZ_HARNESS_H_
