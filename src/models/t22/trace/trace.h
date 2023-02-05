// Copyright 2023 Eliot Courtney.
#ifndef MODELS_T22_TRACE_TRACE_H_
#define MODELS_T22_TRACE_TRACE_H_

#include "compute/energy/t22/model.h"
#include "compute/mfe/t04/dp.h"
#include "compute/traceback/traceback.h"
#include "model/primary.h"

namespace mrna::md::t22::trace {

TraceResult Traceback(
    const Primary& r, const erg::t22::Model::Ptr& em, const DpArray& dp, const ExtArray& ext);

}  // namespace mrna::md::t22::trace

#endif  // MODELS_T22_TRACE_TRACE_H_
