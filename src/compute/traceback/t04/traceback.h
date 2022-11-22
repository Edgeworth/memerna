// Copyright 2021 E.
#ifndef COMPUTE_TRACEBACK_TRACEBACK_H_
#define COMPUTE_TRACEBACK_TRACEBACK_H_

#include <utility>

#include "compute/dp.h"
#include "compute/energy/energy.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::tb::t04 {

TracebackResult Traceback(
    const Primary& r, const energy::EnergyModel& em, const DpArray& dp, const ExtArray& ext);

}  // namespace mrna::tb::t04

#endif  // COMPUTE_TRACEBACK_TRACEBACK_H_
