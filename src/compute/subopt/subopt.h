// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_SUBOPT_SUBOPT_H_
#define COMPUTE_SUBOPT_SUBOPT_H_

#include <functional>

#include "compute/traceback/traceback.h"
#include "model/base.h"
#include "model/ctd.h"

namespace mrna::subopt {

struct SuboptResult {
  SuboptResult() = default;
  SuboptResult(tb::TracebackResult tb, Energy energy) : tb(std::move(tb)), energy(energy) {}

  SuboptResult(SuboptResult&&) = default;
  SuboptResult& operator=(SuboptResult&&) = default;

  // Allow copies explicitly using the constructor.
  explicit SuboptResult(const SuboptResult&) = default;
  SuboptResult& operator=(const SuboptResult&) = delete;

  tb::TracebackResult tb;
  Energy energy;
};

using SuboptCallback = std::function<void(const SuboptResult&)>;

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_SUBOPT_H_
