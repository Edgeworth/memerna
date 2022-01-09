// Copyright 2021 Eliot Courtney.
#ifndef COMPUTE_TRACEBACK_TRACEBACK_H_
#define COMPUTE_TRACEBACK_TRACEBACK_H_

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "model/base.h"

namespace mrna::traceback {

std::tuple<Secondary, Ctds> Traceback(
    const Primary& r, const energy::EnergyModel& em, const DpArray& dp, const ExtArray& ext);

}  // namespace mrna::traceback

#endif  // COMPUTE_TRACEBACK_TRACEBACK_H_
