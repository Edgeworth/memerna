// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_SUBOPT_SUBOPT_H_
#define COMPUTE_SUBOPT_SUBOPT_H_

#include <functional>

#include "model/base.h"
#include "model/globals.h"

namespace mrna::subopt {

typedef std::function<void(const computed_t&)> SuboptimalCallback;

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_SUBOPT_H_
