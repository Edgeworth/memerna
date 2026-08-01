// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BRUTE_ALG_H_
#define BACKENDS_BRUTE_ALG_H_

#include <string>
#include <vector>

#include "api/ctx/backend.h"
#include "api/ctx/backend_cfg.h"
#include "api/energy/energy_cfg.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace_cfg.h"
#include "model/primary.h"

namespace mrna::md::brute {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
class Mfe {
 public:
  [[nodiscard]] static bool IsSupported(
      const trace::TraceCfg& trace_cfg, std::string* reason = nullptr);

  [[nodiscard]] static subopt::SuboptResult Run(const Primary& r, BackendModelPtr m,
      const BackendCfg& backend_cfg, const PseudofreeCfgT& pf, const trace::TraceCfg& trace_cfg);
};

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
class Subopt {
 public:
  [[nodiscard]] static bool IsSupported(
      const subopt::SuboptCfg& subopt_cfg, std::string* reason = nullptr);

  [[nodiscard]] static std::vector<subopt::SuboptResult> Run(const Primary& r, BackendModelPtr m,
      const BackendCfg& backend_cfg, const PseudofreeCfgT& pf, subopt::SuboptCfg subopt_cfg);
};

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
class Pfn {
 public:
  [[nodiscard]] static bool IsSupported(
      const BackendCfg& backend_cfg, std::string* reason = nullptr);

  [[nodiscard]] static pfn::PfnResult Run(
      const Primary& r, BackendModelPtr m, const BackendCfg& backend_cfg, const PseudofreeCfgT& pf);
};

}  // namespace mrna::md::brute

#endif  // BACKENDS_BRUTE_ALG_H_
