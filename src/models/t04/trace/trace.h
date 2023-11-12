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

struct Expansion {
  // Extra energy of this expansion compared to the best choice.
  Energy delta = {ZERO_E};

  // st == -1 used to mean none - using optional here is like a 40% perf hit.
  DpIndex idx0{};
  DpIndex idx1{};
  IndexCtd ctd0{};
  IndexCtd ctd1{};

  bool operator<(const Expansion& o) const { return delta < o.delta; }
};

TraceResult Traceback(
    const Primary& r, const t04::Model::Ptr& em, const trace::TraceCfg& cfg, const DpState& state);

}  // namespace mrna::md::t04

#endif  // MODELS_T04_TRACE_TRACE_H_
