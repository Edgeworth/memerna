// Copyright 2022 E.
#ifndef FUZZ_FUZZ_HARNESS_H_
#define FUZZ_FUZZ_HARNESS_H_

#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "api/energy/model.h"
#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_invocation.h"
#include "model/primary.h"

#ifdef USE_RNASTRUCTURE
#include "api/bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::fuzz {

class FuzzHarness {
 public:
  explicit FuzzHarness(FuzzCfg cfg);

  FuzzInvocation CreateInvocation(
      const Primary& r, std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired);
  std::mt19937& e() { return e_; }

  [[nodiscard]] constexpr uint_fast32_t last_seed() const { return last_seed_; }

 private:
  std::vector<erg::EnergyModelPtr> ems_;
  uint_fast32_t last_seed_ = 0;
  std::vector<Energy> last_pf_paired_;
  std::vector<Energy> last_pf_unpaired_;
  FuzzCfg cfg_;
  std::mt19937 e_;

#ifdef USE_RNASTRUCTURE
  std::shared_ptr<bridge::RNAstructure> rstr_;
#endif  // USE_RNASTRUCTURE

  void MaybeLoadModels(
      uint_fast32_t seed, std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired);
  [[nodiscard]] erg::EnergyModelPtr LoadModel(
      const std::string& model_name, uint_fast32_t seed) const;
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZ_HARNESS_H_
