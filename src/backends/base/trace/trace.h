// Copyright 2021 Eliot Courtney.
#ifndef BACKENDS_BASE_TRACE_TRACE_H_
#define BACKENDS_BASE_TRACE_TRACE_H_

#include <string>

#include "api/energy/energy_cfg.h"
#include "api/trace/trace.h"
#include "api/trace/trace_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
class Trace {
 public:
  static bool IsSupported(
      const PseudofreeCfgT& pf, const trace::TraceCfg& tcfg, std::string* reason = nullptr);

  static trace::TraceResult Run(const Primary& r, const Model::Ptr& m, const DpState& state,
      const PseudofreeCfgT& pf, const trace::TraceCfg& tcfg);
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_TRACE_TRACE_H_
