// Copyright 2022 E.
#include "compute/brute/alg.h"

#include <set>
#include <utility>

#include "compute/brute/brute.h"
#include "compute/subopt/config.h"

namespace mrna::brute {

subopt::SuboptResult MfeBruteForce(const Primary& r, energy::EnergyModelPtr em) {
  return std::move(SuboptimalBruteForce(r, em, subopt::SuboptCfg{.strucs = 1})[0]);
}

part::PartResult PartitionBruteForce(const Primary& r, energy::EnergyModelPtr em) {
  // TODO: Allow lonely pairs for the partition function
  auto res = BruteForce(r, em, {.part = true}).Run();
  return {.dp{}, .ext{}, .part = std::move(res.part), .prob = std::move(res.prob)};
}

std::vector<subopt::SuboptResult> SuboptimalBruteForce(
    const Primary& r, energy::EnergyModelPtr em, subopt::SuboptCfg cfg) {
  auto res = BruteForce(r, em, {.subopt = true, .subopt_cfg = cfg}).Run();
  return std::vector<subopt::SuboptResult>(res.subopts.begin(), res.subopts.end());
}

}  // namespace mrna::brute
