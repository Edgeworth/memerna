// Copyright 2023 E.
#ifndef MODELS_T22_TRACE_TRACE_H_
#define MODELS_T22_TRACE_TRACE_H_

#include "api/trace.h"
#include "model/primary.h"
#include "models/t22/energy/model.h"
#include "models/t22/mfe/mfe.h"

namespace mrna::md::t22 {

using trace::TraceResult;

TraceResult Traceback(const Primary& r, const Model::Ptr& em, const DpState& state);

}  // namespace mrna::md::t22

#endif  // MODELS_T22_TRACE_TRACE_H_
