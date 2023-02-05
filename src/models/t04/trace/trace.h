// Copyright 2021 Eliot Courtney.
#ifndef COMPUTE_TRACEBACK_T04_TRACEBACK_H_
#define COMPUTE_TRACEBACK_T04_TRACEBACK_H_

#include "compute/traceback/traceback.h"
#include "model/primary.h"
#include "models/t04/energy/model.h"
#include "models/t04/mfe/dp.h"

namespace mrna::md::t04::trace {

TraceResult Traceback(
    const Primary& r, const erg::t04::Model::Ptr& em, const DpArray& dp, const ExtArray& ext);

}  // namespace mrna::md::t04::trace

#endif  // COMPUTE_TRACEBACK_T04_TRACEBACK_H_
