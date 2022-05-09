// Copyright 2022 E.
#ifndef FUZZ_FUZZ_INVOCATION_H_
#define FUZZ_FUZZ_INVOCATION_H_
#include <deque>
#include <memory>
#include <string>
#include <vector>

#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/config.h"
#include "compute/subopt/subopt.h"
#include "ctx/ctx.h"
#include "fuzz/config.h"
#include "model/primary.h"

#ifdef USE_RNASTRUCTURE
#include <utility>

#include "bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::fuzz {

using Error = std::deque<std::string>;

class FuzzInvocation {
 public:
  FuzzInvocation(const Primary& r, energy::EnergyModelPtr em, const FuzzCfg& cfg);

  Error Run();

#ifdef USE_RNASTRUCTURE
  void set_rnastructure(std::shared_ptr<bridge::RNAstructure> rstr) { rstr_ = std::move(rstr); }
#endif  // USE_RNASTRUCTURE

 private:
  Primary r_;
  energy::EnergyModelPtr em_;
  FuzzCfg cfg_;

  ctx::FoldResult fold_{};
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