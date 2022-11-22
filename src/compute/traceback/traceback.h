// Copyright 2021 Eliot Courtney.
#ifndef COMPUTE_TRACEBACK_TRACEBACK_H_
#define COMPUTE_TRACEBACK_TRACEBACK_H_

#include <utility>

#include "compute/dp.h"
#include "compute/energy/energy.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::tb {

struct TracebackResult {
  TracebackResult() = default;
  ~TracebackResult() = default;
  TracebackResult(Secondary s, Ctds ctd) : s(std::move(s)), ctd(std::move(ctd)) {}

  TracebackResult(TracebackResult&&) = default;
  TracebackResult& operator=(TracebackResult&&) = default;

  // Allow copies explicitly using the constructor.
  explicit TracebackResult(const TracebackResult&) = default;
  TracebackResult& operator=(const TracebackResult&) = delete;

  constexpr auto operator<=>(const TracebackResult&) const = default;

  Secondary s;
  Ctds ctd;  // May be empty if CTDs were not computed.
};

}  // namespace mrna::tb

#endif  // COMPUTE_TRACEBACK_TRACEBACK_H_
