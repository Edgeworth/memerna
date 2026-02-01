// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BRUTE_ALG_H_
#define BACKENDS_BRUTE_ALG_H_

#include <vector>

#include "api/ctx/backend.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "model/primary.h"

namespace mrna::md::brute {

subopt::SuboptResult MfeBrute(
    const Primary& r, BackendModelPtr m, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf);

pfn::PfnResult PfnBrute(
    const Primary& r, BackendModelPtr m, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf);

std::vector<subopt::SuboptResult> SuboptBrute(const Primary& r, BackendModelPtr m,
    erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, subopt::SuboptCfg subopt_cfg);

}  // namespace mrna::md::brute

#endif  // BACKENDS_BRUTE_ALG_H_
