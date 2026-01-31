// Copyright 2021 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_TRACE_TRACE_H_
#define BACKENDS_BASEOPT_TRACE_TRACE_H_

#include "api/energy/pseudofree_cfg.h"
#include "api/trace/trace_cfg.h"
#include "backends/baseopt/energy/model.h"
#include "backends/common/base/dp.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

TraceResult Traceback(const Primary& r, const Model::Ptr& m, const DpState& state,
    const erg::PseudofreeCfg& pf, const trace::TraceCfg& cfg);

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_TRACE_TRACE_H_
