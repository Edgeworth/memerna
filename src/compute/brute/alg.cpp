// Copyright 2022 Eliot Courtney.
#include "compute/brute/alg.h"

#include <set>
#include <utility>

#include "compute/brute/brute.h"
#include "compute/subopt/config.h"

namespace mrna::brute {

subopt::SuboptResult MfeBruteForce(const Primary& r, energy::EnergyModelPtr em) {
  return std::move(SuboptimalBruteForce(r, std::move(em), subopt::SuboptCfg{.strucs = 1})[0]);
}

part::PartResult PartitionBruteForce(const Primary& r, energy::EnergyModelPtr em) {
  // TODO: Allow lonely pairs for the partition function
  auto res = BruteForce(r, std::move(em), {.part = true}).Run();
  return {.dp{}, .ext{}, .part = std::move(res.part), .prob = std::move(res.prob)};
}

std::vector<subopt::SuboptResult> SuboptimalBruteForce(
    const Primary& r, energy::EnergyModelPtr em, subopt::SuboptCfg cfg) {
  auto res = BruteForce(r, std::move(em), {.subopt = true, .subopt_cfg = cfg}).Run();
  return {res.subopts.begin(), res.subopts.end()};
}

}  // namespace mrna::brute
