// Copyright 2021 E.
#ifndef COMPUTE_TRACEBACK_TRACEBACK_H_
#define COMPUTE_TRACEBACK_TRACEBACK_H_

#include <optional>

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "model/base.h"

namespace mrna::traceback {

struct TracebackResult {
  Secondary s;
  Ctds ctd;  // May be empty if CTDs were not computed.
};

TracebackResult Traceback(
    const Primary& r, const energy::EnergyModel& em, const DpArray& dp, const ExtArray& ext);

}  // namespace mrna::traceback

#endif  // COMPUTE_TRACEBACK_TRACEBACK_H_
