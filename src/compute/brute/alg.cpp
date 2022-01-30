// Copyright 2022 E.
#include "compute/brute/alg.h"

#include <set>
#include <utility>

#include "compute/brute/brute.h"
#include "compute/subopt/config.h"

namespace mrna::brute {

subopt::SuboptResult MfeBruteForce(Primary r, const energy::EnergyModel& em) {
  return std::move(SuboptimalBruteForce(std::move(r), em, subopt::SuboptCfg{.strucs = 1})[0]);
}

part::PartResult PartitionBruteForce(Primary r, const energy::EnergyModel& em) {
  // TODO: Allow lonely pairs for the partition function
  auto res = BruteForce(std::move(r), em, {.part = true}).Run();
  return {.dp{}, .ext{}, .part = std::move(res.part), .prob = std::move(res.prob)};
}

std::vector<subopt::SuboptResult> SuboptimalBruteForce(
    Primary r, const energy::EnergyModel& em, subopt::SuboptCfg cfg) {
  auto res = BruteForce(std::move(r), em, {.subopt = true, .subopt_cfg = cfg}).Run();
  return std::vector<subopt::SuboptResult>(res.subopts.begin(), res.subopts.end());
}

}  // namespace mrna::brute
