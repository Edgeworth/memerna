// Copyright 2022 Eliot Courtney.
#include "backends/brute/alg.h"

#include <set>
#include <utility>
#include <vector>

#include "api/subopt/subopt_cfg.h"
#include "backends/brute/brute.h"
#include "backends/common/static.h"

namespace mrna::md::brute {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
subopt::SuboptResult MfeBrute(
    const Primary& r, BackendModelPtr m, const BackendCfg& backend_cfg, const PseudofreeCfgT& pf) {
  return std::move(SuboptBrute<CFG, PseudofreeCfgT>(
      r, std::move(m), backend_cfg, pf, subopt::SuboptCfg{.strucs = 1})[0]);
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
pfn::PfnResult PfnBrute(
    const Primary& r, BackendModelPtr m, const BackendCfg& backend_cfg, const PseudofreeCfgT& pf) {
  // TODO(2): Allow lonely pairs for the partition function
  auto res = Brute<CFG, PseudofreeCfgT>(r, std::move(m), backend_cfg, pf, {.pfn = true}).Run();
  return {.state{}, .pfn = std::move(res.pfn)};
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
std::vector<subopt::SuboptResult> SuboptBrute(const Primary& r, BackendModelPtr m,
    const BackendCfg& backend_cfg, const PseudofreeCfgT& pf, subopt::SuboptCfg subopt_cfg) {
  auto res = Brute<CFG, PseudofreeCfgT>(
      r, std::move(m), backend_cfg, pf, {.subopt = true, .subopt_cfg = subopt_cfg})
                 .Run();
  return {res.subopts.begin(), res.subopts.end()};
}

#define MRNA_INSTANTIATE(CFG, PF)                                     \
  template subopt::SuboptResult MfeBrute<CFG, PF>(                    \
      const Primary&, BackendModelPtr, const BackendCfg&, const PF&); \
  template pfn::PfnResult PfnBrute<CFG, PF>(                          \
      const Primary&, BackendModelPtr, const BackendCfg&, const PF&); \
  template std::vector<subopt::SuboptResult> SuboptBrute<CFG, PF>(    \
      const Primary&, BackendModelPtr, const BackendCfg&, const PF&, subopt::SuboptCfg);
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::brute
