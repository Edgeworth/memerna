// Copyright 2021 E.
#ifndef COMPUTE_TRACEBACK_T04_TRACEBACK_H_
#define COMPUTE_TRACEBACK_T04_TRACEBACK_H_

#include "api/trace.h"
#include "model/primary.h"
#include "models/t04/energy/model.h"
#include "models/t04/mfe/dp.h"

namespace mrna::md::t04 {

using trace::TraceResult;

TraceResult Traceback(const Primary& r, const t04::Model::Ptr& em, const DpState& state);

}  // namespace mrna::md::t04

#endif  // COMPUTE_TRACEBACK_T04_TRACEBACK_H_
