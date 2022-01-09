// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_SUBOPT_SUBOPT_H_
#define COMPUTE_SUBOPT_SUBOPT_H_

#include <functional>

#include "model/base.h"
#include "model/ctd.h"

namespace mrna::subopt {

using SuboptimalCallback = std::function<void(const Computed&)>;

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_SUBOPT_H_
