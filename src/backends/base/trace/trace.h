// Copyright 2021 Eliot Courtney.
#ifndef BACKENDS_BASE_TRACE_TRACE_H_
#define BACKENDS_BASE_TRACE_TRACE_H_

#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/trace/trace_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

TraceResult Traceback(const Primary& r, const base::Model::Ptr& m, const DpState& state,
    erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, const trace::TraceCfg& tcfg);

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_TRACE_TRACE_H_
