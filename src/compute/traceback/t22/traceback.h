// Copyright 2023 Eliot Courtney.
#ifndef COMPUTE_TRACEBACK_T22_TRACEBACK_H_
#define COMPUTE_TRACEBACK_T22_TRACEBACK_H_

#include "compute/dp.h"
#include "compute/energy/t22/model.h"
#include "compute/traceback/traceback.h"
#include "model/primary.h"

namespace mrna::tb::t22 {

TracebackResult Traceback(
    const Primary& r, const erg::t22::Model::Ptr& em, const DpArray& dp, const ExtArray& ext);

}  // namespace mrna::tb::t22

#endif  // COMPUTE_TRACEBACK_T22_TRACEBACK_H_
