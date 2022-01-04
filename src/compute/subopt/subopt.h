// Copyright 2016 E.
#ifndef COMPUTE_SUBOPT_SUBOPT_H_
#define COMPUTE_SUBOPT_SUBOPT_H_

#include <functional>

#include "model/base.h"
#include "model/ctd.h"
#include "model/globals.h"

namespace mrna::subopt {

typedef std::function<void(const Computed&)> SuboptimalCallback;

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_SUBOPT_H_
