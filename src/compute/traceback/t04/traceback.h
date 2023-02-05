// Copyright 2021 Eliot Courtney.
#ifndef COMPUTE_TRACEBACK_T04_TRACEBACK_H_
#define COMPUTE_TRACEBACK_T04_TRACEBACK_H_

#include "compute/energy/t04/model.h"
#include "compute/mfe/t04/dp.h"
#include "compute/traceback/traceback.h"
#include "model/primary.h"

namespace mrna::tb::t04 {

TracebackResult Traceback(
    const Primary& r, const erg::t04::Model::Ptr& em, const DpArray& dp, const ExtArray& ext);

}  // namespace mrna::tb::t04

#endif  // COMPUTE_TRACEBACK_T04_TRACEBACK_H_
