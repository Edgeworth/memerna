// Copyright 2021 E.
#ifndef FUZZ_FUZZER_H_
#define FUZZ_FUZZER_H_

#include <deque>
#include <string>
#include <vector>

#include "compute/dp.h"
#include "compute/energy/model.h"
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
  Fuzzer(Primary r, energy::EnergyModel em, FuzzCfg cfg);

  Error Run();

 private:
  Primary r_;
  energy::EnergyModel em_;
  FuzzCfg cfg_;

  Error errors_;

  std::vector<subopt::SuboptResult> memerna_subopts_;
  std::vector<DpArray> memerna_dps;

  // RNAstructure related:
#ifdef USE_RNASTRUCTURE
  std::optional<bridge::RNAstructure> rnastructure_;
  dp_state_t rnastructure_dp_;

  void set_rnastructure(bridge::RNAstructure&& rnastructure) {
    rnastructure_ = std::move(rnastructure);
  }
#endif  // USE_RNASTRUCTURE

  void Register(const std::string& header, Error&& local);

  bool HasDuplicates(const std::vector<subopt::SuboptResult>& subopts);

  Error CheckSuboptimalResult(const std::vector<subopt::SuboptResult>& subopt, bool has_ctds);

  Error CheckSuboptimalResultPair(
      const std::vector<subopt::SuboptResult>& a, const std::vector<subopt::SuboptResult>& b);

  Error CheckSuboptimal();

  Error CheckDpTables();

  Error MemernaComputeAndCheckState();

  Error RnastructureComputeAndCheckState();

  Error CheckBruteForce();

  Error CheckPartition();
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZER_H_
