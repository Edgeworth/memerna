// Copyright 2021 Eliot Courtney.
#ifndef FUZZ_FUZZER_H_
#define FUZZ_FUZZER_H_

#include <deque>
#include <optional>

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "fuzz/config.h"
#include "model/structure.h"
#include "util/array.h"

#ifdef USE_RNASTRUCTURE
#include "bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::fuzz {

typedef std::deque<std::string> Error;

class Fuzzer {
 public:
  Fuzzer(Primary r_, const FuzzCfg& cfg_, const energy::EnergyModelPtr em_);

  Error Run();

 private:
  const int N;
  const Primary r_;
  const FuzzCfg cfg_;
  const energy::EnergyModelPtr em_;

  std::vector<Computed> memerna_computeds_;
  std::vector<Array3D<Energy, DP_SIZE>> memerna_dps;

  // RNAstructure related:
#ifdef USE_RNASTRUCTURE
  std::optional<bridge::RNAstructure> rnastructure_;
  dp_state_t rnastructure_dp_;

  void set_rnastructure(bridge::RNAstructure&& rnastructure) {
    rnastructure_ = std::move(rnastructure);
  }
#endif  // USE_RNASTRUCTURE

  Error MaybePrepend(const Error& main, const std::string& header);

  void AppendErrors(Error& main, Error&& extra);

  bool HasDuplicates(const std::vector<Computed>& computeds);

  Error CheckSuboptimalResult(const std::vector<Computed>& subopt, bool has_ctds);

  Error CheckSuboptimalResultPair(const std::vector<Computed>& a, const std::vector<Computed>& b);

  Error CheckSuboptimal();

  Error CheckDpTables();

  Error MemernaComputeAndCheckState();

  Error RnastructureComputeAndCheckState();

  Error CheckBruteForce();

  Error CheckPartition();
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZER_H_
