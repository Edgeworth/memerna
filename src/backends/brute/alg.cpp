// Copyright 2022 Eliot Courtney.
#include "backends/brute/alg.h"

#include <set>
#include <utility>
#include <vector>

#include "api/subopt/subopt_cfg.h"
#include "backends/brute/brute.h"
#include "backends/common/static.h"
#include "util/error.h"

namespace mrna::md::brute {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
bool Mfe<CFG, PseudofreeCfgT>::IsSupported(const trace::TraceCfg& trace_cfg, std::string* reason) {
  if (!trace_cfg.random) return true;
  if (reason) *reason = "random traceback is not supported by MFE brute";
  return false;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
bool Subopt<CFG, PseudofreeCfgT>::IsSupported(
    const subopt::SuboptCfg& subopt_cfg, std::string* reason) {
  if (subopt_cfg.time_secs < 0.0F) return true;
  if (reason) *reason = "brute force does not support a time limit";
  return false;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
bool Pfn<CFG, PseudofreeCfgT>::IsSupported(const BackendCfg& backend_cfg, std::string* reason) {
  if (backend_cfg.energy_model != erg::EnergyModelKind::T22) return true;
  if (reason) *reason = "brute partition function does not support T22";
  return false;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
subopt::SuboptResult Mfe<CFG, PseudofreeCfgT>::Run(const Primary& r, BackendModelPtr m,
    const BackendCfg& backend_cfg, const PseudofreeCfgT& pf, const trace::TraceCfg& trace_cfg) {
  trace_cfg.Verify();
  std::string reason;
  verify(IsSupported(trace_cfg, &reason), "MFE brute does not support the given configuration: {}",
      reason);
  return std::move(Subopt<CFG, PseudofreeCfgT>::Run(
      r, std::move(m), backend_cfg, pf, subopt::SuboptCfg{.strucs = 1})[0]);
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
pfn::PfnResult Pfn<CFG, PseudofreeCfgT>::Run(
    const Primary& r, BackendModelPtr m, const BackendCfg& backend_cfg, const PseudofreeCfgT& pf) {
  std::string reason;
  verify(IsSupported(backend_cfg, &reason),
      "PFN brute does not support the given configuration: {}", reason);
  // TODO(2): Allow lonely pairs for the partition function
  auto res = Brute<CFG, PseudofreeCfgT>(r, std::move(m), backend_cfg, pf, {.pfn = true}).Run();
  return {.state{}, .pfn = std::move(res.pfn)};
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
std::vector<subopt::SuboptResult> Subopt<CFG, PseudofreeCfgT>::Run(const Primary& r,
    BackendModelPtr m, const BackendCfg& backend_cfg, const PseudofreeCfgT& pf,
    subopt::SuboptCfg subopt_cfg) {
  subopt_cfg.Verify();
  std::string reason;
  verify(IsSupported(subopt_cfg, &reason),
      "subopt brute does not support the given configuration: {}", reason);
  auto res = Brute<CFG, PseudofreeCfgT>(
      r, std::move(m), backend_cfg, pf, {.subopt = true, .subopt_cfg = subopt_cfg})
                 .Run();
  return {res.subopts.begin(), res.subopts.end()};
}

#define MRNA_INSTANTIATE(CFG, PF) \
  template class Mfe<CFG, PF>;    \
  template class Subopt<CFG, PF>; \
  template class Pfn<CFG, PF>;
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::brute
