// Copyright 2021 Eliot Courtney.
#ifndef BACKENDS_BASE_TRACE_TRACE_H_
#define BACKENDS_BASE_TRACE_TRACE_H_

#include "api/trace/trace_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

TraceResult Traceback(
    const Primary& r, const base::Model::Ptr& m, const trace::TraceCfg& cfg, const DpState& state);

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_TRACE_TRACE_H_
