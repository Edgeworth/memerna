// Copyright 2023 Eliot Courtney.
#ifndef BACKENDS_STACK_TRACE_TRACE_H_
#define BACKENDS_STACK_TRACE_TRACE_H_

#include "api/trace/trace.h"
#include "api/trace/trace_cfg.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/mfe.h"
#include "model/primary.h"

namespace mrna::md::stack {

using trace::TraceResult;

// Holds information about a particular recurrence (expansion).
struct Expansion {
  // Extra energy of this expansion compared to the best choice (not used for
  // traceback, it will always be zero - used for suboptimal folding).
  Energy delta = {ZERO_E};

  std::optional<DpIndex> idx0 = std::nullopt;
  std::optional<DpIndex> idx1 = std::nullopt;
  IndexCtd ctd0{};
  IndexCtd ctd1{};
  Pair pair{};

  bool operator<(const Expansion& o) const { return delta < o.delta; }
};

TraceResult Traceback(
    const Primary& r, const Model::Ptr& m, const trace::TraceCfg& cfg, const DpState& state);

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_TRACE_TRACE_H_
