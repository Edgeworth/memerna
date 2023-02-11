// Copyright 2022 Eliot Courtney.
#ifndef FUZZ_FUZZ_INVOCATION_H_
#define FUZZ_FUZZ_INVOCATION_H_
#include <deque>
#include <string>
#include <vector>

#include "api/ctx/ctx.h"
#include "api/energy/model.h"
#include "api/part.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "fuzz/fuzz_cfg.h"
#include "model/primary.h"

#ifdef USE_RNASTRUCTURE
#include <memory>
#include <utility>

#include "api/bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::fuzz {

using Error = std::deque<std::string>;

class FuzzInvocation {
 public:
  FuzzInvocation(const Primary& r, erg::EnergyModelPtr em, const FuzzCfg& cfg, uint_fast32_t seed);

  Error Run();

  [[nodiscard]] uint_fast32_t seed() const { return seed_; }

#ifdef USE_RNASTRUCTURE
  void set_rnastructure(std::shared_ptr<bridge::RNAstructure> rstr) { rstr_ = std::move(rstr); }
#endif  // USE_RNASTRUCTURE

 private:
  Primary r_;
  erg::EnergyModelPtr em_;
  FuzzCfg cfg_;
  uint_fast32_t seed_;

  std::optional<FoldResult> fold_;
  std::vector<subopt::SuboptResult> subopt_{};
  part::PartResult part_{};
  Error errors_;

#ifdef USE_RNASTRUCTURE
  std::shared_ptr<bridge::RNAstructure> rstr_;

  Error CheckMfeRNAstructure();
  Error CheckSuboptRNAstructure(subopt::SuboptCfg cfg);
  Error CheckPartitionRNAstructure();
#endif  // USE_RNASTRUCTURE

  void Register(const std::string& header, Error&& local);

  void EnsureFoldResult();

  Error CheckMfe();

  Error CheckSubopt();

  static bool SuboptDuplicates(const std::vector<subopt::SuboptResult>& subopts);
  Error CheckSuboptResult(const std::vector<subopt::SuboptResult>& subopt, bool has_ctds);
  static Error CheckSuboptResultPair(subopt::SuboptCfg cfg,
      const std::vector<subopt::SuboptResult>& a, const std::vector<subopt::SuboptResult>& b);

  Error CheckPartition();
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZ_INVOCATION_H_
