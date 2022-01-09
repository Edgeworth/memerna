// Copyright 2016 E.
#ifndef COMPUTE_SUBOPT_SUBOPT_H_
#define COMPUTE_SUBOPT_SUBOPT_H_

#include <functional>

#include "compute/traceback/traceback.h"
#include "model/base.h"
#include "model/ctd.h"

namespace mrna::subopt {

struct SuboptResult {
  traceback::TracebackResult tb;
  Energy energy;
};

using SuboptCallback = std::function<void(const SuboptResult&)>;

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_SUBOPT_H_
