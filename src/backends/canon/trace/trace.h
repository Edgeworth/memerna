// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_TRACE_TRACE_H_
#define BACKENDS_CANON_TRACE_TRACE_H_

#include <string>

#include "api/energy/energy_cfg.h"
#include "api/trace/trace.h"
#include "api/trace/trace_cfg.h"
#include "backends/canon/dp.h"
#include "backends/canon/energy/model.h"
#include "model/primary.h"

namespace mrna::md::canon {

using trace::TraceResult;

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
class Trace {
 public:
  static bool IsSupported(
      const PseudofreeCfgT& pf, const trace::TraceCfg& tcfg, std::string* reason = nullptr);

  static TraceResult Run(const Primary& r, const Model::Ptr& m, const DpState& state,
      const PseudofreeCfgT& pf, const trace::TraceCfg& tcfg);
};

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_TRACE_TRACE_H_
