// Copyright 2016 Eliot Courtney.
#include "api/bridge/memerna.h"

#include <cstdint>
#include <optional>
#include <vector>

#include "api/ctx/algorithm.h"
#include "api/ctx/ctx.h"
#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace_cfg.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::bridge {

erg::EnergyResult Memerna::Efn(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, const Ctds* given_ctd, bool build_structure) const {
  return ctx_.Efn(r, s, cfg, pf, given_ctd, build_structure);
}

FoldResult Memerna::Fold(const Primary& r, std::optional<MfeAlg> alg, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) const {
  return ctx_.Fold(r, alg, cfg, pf, trace_cfg);
}

int64_t Memerna::Subopt(const Primary& r, std::optional<MfeAlg> mfe_alg,
    std::optional<SuboptAlg> alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
    subopt::SuboptCfg subopt_cfg, const subopt::SuboptCallback& fn) const {
  return ctx_.Subopt(r, mfe_alg, alg, cfg, pf, subopt_cfg, fn);
}

std::vector<subopt::SuboptResult> Memerna::SuboptIntoVector(const Primary& r,
    std::optional<MfeAlg> mfe_alg, std::optional<SuboptAlg> alg, erg::EnergyCfg cfg,
    const erg::PseudofreeCfg& pf, subopt::SuboptCfg subopt_cfg) const {
  return ctx_.SuboptIntoVector(r, mfe_alg, alg, cfg, pf, subopt_cfg);
}

pfn::PfnResult Memerna::Pfn(const Primary& r, std::optional<PfnAlg> alg, erg::EnergyCfgPfn cfg,
    const erg::PseudofreeCfg& pf) const {
  return ctx_.Pfn(r, alg, cfg, pf);
}

Memerna Memerna::FromArgParse(const ArgParse& args) { return Memerna(Ctx::FromArgParse(args)); }

}  // namespace mrna::bridge
