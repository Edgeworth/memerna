// Copyright 2021 Eliot Courtney.
#ifndef FUZZ_FUZZER_H_
#define FUZZ_FUZZER_H_

#include <deque>
#include <memory>
#include <string>
#include <vector>

#include "bridge/bridge.h"
#include "compute/dp.h"
#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "fuzz/config.h"
#include "model/primary.h"

#ifdef USE_RNASTRUCTURE
#include <utility>

#include "bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::fuzz {

using Error = std::deque<std::string>;

class Fuzzer {
 public:
  Fuzzer(Primary r, energy::EnergyModelPtr em, FuzzCfg cfg);

  Error Run();

 private:
  Primary r_;
  energy::EnergyModelPtr em_;
  FuzzCfg cfg_;

  ctx::FoldResult fold_;
  std::vector<subopt::SuboptResult> subopt_;
  part::PartResult part_;
  Error errors_;

#ifdef USE_RNASTRUCTURE
  std::shared_ptr<bridge::RNAstructure> rstr_;
  void set_rnastructure(std::shared_ptr<bridge::RNAstructure> rstr) { rstr_ = std::move(rstr); }
#endif  // USE_RNASTRUCTURE

  void Register(const std::string& header, Error&& local);

  Error CheckMfe();
  Error CheckMfeRNAstructure();

  Error CheckSubopt();
  Error CheckSuboptRNAstructure();

  bool SuboptDuplicates(const std::vector<subopt::SuboptResult>& subopts);
  Error CheckSuboptResult(const std::vector<subopt::SuboptResult>& subopt, bool has_ctds);
  Error CheckSuboptResultPair(
      const std::vector<subopt::SuboptResult>& a, const std::vector<subopt::SuboptResult>& b);

  Error CheckPartition();
  Error CheckPartitionRNAstructure();
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZER_H_
