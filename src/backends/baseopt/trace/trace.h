// Copyright 2021 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_TRACE_TRACE_H_
#define BACKENDS_BASEOPT_TRACE_TRACE_H_

#include "api/trace/trace_cfg.h"
#include "backends/baseopt/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

TraceResult Traceback(
    const Primary& r, const Model::Ptr& m, const trace::TraceCfg& cfg, const DpState& state);

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_TRACE_TRACE_H_
