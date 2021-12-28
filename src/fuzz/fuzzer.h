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

typedef std::deque<std::string> error_t;

class Fuzzer {
 public:
  Fuzzer(primary_t r_, const cfg_t& cfg_, const energy::EnergyModelPtr em_);

  error_t Run();

 private:
  const int N;
  const primary_t r;
  const cfg_t cfg;
  const energy::EnergyModelPtr em;

  std::vector<computed_t> memerna_computeds;
  std::vector<array3d_t<energy_t, DP_SIZE>> memerna_dps;

  // RNAstructure related:
#ifdef USE_RNASTRUCTURE
  std::optional<bridge::RNAstructure> rnastructure;
  dp_state_t rnastructure_dp;

  void set_rnastructure(bridge::RNAstructure&& rnastructure) {
    rnastructure = std::move(rnastructure);
  }
#endif  // USE_RNASTRUCTURE

  error_t MaybePrepend(const error_t& main, const std::string& header);

  void AppendErrors(error_t& main, error_t&& extra);

  bool HasDuplicates(const std::vector<computed_t>& computeds);

  error_t CheckSuboptimalResult(const std::vector<computed_t>& subopt, bool has_ctds);

  error_t CheckSuboptimalResultPair(
      const std::vector<computed_t>& a, const std::vector<computed_t>& b);

  error_t CheckSuboptimal();

  error_t CheckDpTables();

  error_t MemernaComputeAndCheckState();

  error_t RnastructureComputeAndCheckState();

  error_t CheckBruteForce();

  error_t CheckPartition();
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZER_H_
