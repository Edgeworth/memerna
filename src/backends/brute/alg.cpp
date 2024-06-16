// Copyright 2022 Eliot Courtney.
#include "backends/brute/alg.h"

#include <set>
#include <utility>

#include "api/subopt/subopt_cfg.h"
#include "backends/brute/brute.h"

namespace mrna::md::brute {

subopt::SuboptResult MfeBrute(const Primary& r, BackendModelPtr m) {
  return std::move(SuboptBrute(r, std::move(m), subopt::SuboptCfg{.strucs = 1})[0]);
}

pfn::PfnResult PfnBrute(const Primary& r, BackendModelPtr m) {
  // TODO(2): Allow lonely pairs for the partition function
  auto res = Brute(r, std::move(m), {.pfn = true}).Run();
  return {.state{}, .pfn = std::move(res.pfn)};
}

std::vector<subopt::SuboptResult> SuboptBrute(
    const Primary& r, BackendModelPtr m, subopt::SuboptCfg cfg) {
  auto res = Brute(r, std::move(m), {.subopt = true, .subopt_cfg = cfg}).Run();
  return {res.subopts.begin(), res.subopts.end()};
}

}  // namespace mrna::md::brute
