// Copyright 2022 Eliot Courtney.
#ifndef FUZZ_FUZZ_HARNESS_H_
#define FUZZ_FUZZ_HARNESS_H_

#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_invocation.h"
#include "model/primary.h"

#ifdef USE_RNASTRUCTURE
#include "api/bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::fuzz {

class FuzzHarness {
 public:
  explicit FuzzHarness(FuzzCfg fuzz_cfg);

  FuzzInvocation CreateInvocation(
      const Primary& r, std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired);
  std::mt19937& e() { return e_; }

  [[nodiscard]] constexpr std::optional<uint_fast32_t> last_seed() const {
    return backend_cfg_.seed;
  }

 private:
  std::vector<BackendModelPtr> ms_;
  BackendCfg backend_cfg_{};
  FuzzCfg fuzz_cfg_;
  std::mt19937 e_;

#ifdef USE_RNASTRUCTURE
  std::shared_ptr<bridge::RNAstructure> rstr_;
#endif  // USE_RNASTRUCTURE

  void MaybeLoadBackends(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired);
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZ_HARNESS_H_
