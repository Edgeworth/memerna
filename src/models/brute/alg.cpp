// Copyright 2022 Eliot Courtney.
#include "models/brute/alg.h"

#include <set>
#include <utility>

#include "api/subopt/subopt_cfg.h"
#include "models/brute/brute.h"

namespace mrna::md::brute {

subopt::SuboptResult MfeBrute(const Primary& r, erg::EnergyModelPtr em) {
  return std::move(SuboptBrute(r, std::move(em), subopt::SuboptCfg{.strucs = 1})[0]);
}

part::PartResult PartitionBrute(const Primary& r, erg::EnergyModelPtr em) {
  // TODO(2): Allow lonely pairs for the partition function
  auto res = Brute(r, std::move(em), {.part = true}).Run();
  return {.state{}, .part = std::move(res.part)};
}

std::vector<subopt::SuboptResult> SuboptBrute(
    const Primary& r, erg::EnergyModelPtr em, subopt::SuboptCfg cfg) {
  auto res = Brute(r, std::move(em), {.subopt = true, .subopt_cfg = cfg}).Run();
  return {res.subopts.begin(), res.subopts.end()};
}

}  // namespace mrna::md::brute
