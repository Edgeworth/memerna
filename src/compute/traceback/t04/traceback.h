// Copyright 2021 E.
#ifndef COMPUTE_TRACEBACK_T04_TRACEBACK_H_
#define COMPUTE_TRACEBACK_T04_TRACEBACK_H_

#include "compute/dp.h"
#include "compute/energy/t04/model.h"
#include "compute/traceback/traceback.h"
#include "model/primary.h"

namespace mrna::tb::t04 {

TracebackResult Traceback(
    const Primary& r, const erg::t04::ModelPtr& em, const DpArray& dp, const ExtArray& ext);

}  // namespace mrna::tb::t04

#endif  // COMPUTE_TRACEBACK_T04_TRACEBACK_H_