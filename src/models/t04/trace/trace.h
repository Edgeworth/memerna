// Copyright 2021 Eliot Courtney.
#ifndef MODELS_T04_TRACE_TRACE_H_
#define MODELS_T04_TRACE_TRACE_H_

#include "api/trace/trace.h"
#include "api/trace/trace_cfg.h"
#include "model/primary.h"
#include "models/t04/energy/model.h"
#include "models/t04/mfe/dp.h"

namespace mrna::md::t04 {

using trace::TraceResult;

TraceResult Traceback(
    const Primary& r, const t04::Model::Ptr& em, const trace::TraceCfg& cfg, const DpState& state);

}  // namespace mrna::md::t04

#endif  // MODELS_T04_TRACE_TRACE_H_
