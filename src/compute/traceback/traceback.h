// Copyright 2021 E.
#ifndef COMPUTE_TRACEBACK_TRACEBACK_H_
#define COMPUTE_TRACEBACK_TRACEBACK_H_

#include "compute/dp.h"
#include "compute/energy/model.h"
#include "model/base.h"

namespace mrna::traceback {

std::tuple<std::vector<int>, std::vector<Ctd>> Traceback(
    const Primary& r, const energy::EnergyModel& em, const DpArray& dp, const ExtArray& ext);

}  // namespace mrna::traceback

#endif  // COMPUTE_TRACEBACK_TRACEBACK_H_
