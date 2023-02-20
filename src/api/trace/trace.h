// Copyright 2021 Eliot Courtney.
#ifndef API_TRACE_TRACE_H_
#define API_TRACE_TRACE_H_

#include <utility>

#include "api/energy/model.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/t04/mfe/dp.h"

namespace mrna::trace {

struct TraceResult {
  TraceResult() = default;
  ~TraceResult() = default;
  TraceResult(Secondary s, Ctds ctd) : s(std::move(s)), ctd(std::move(ctd)) {}

  TraceResult(TraceResult&&) = default;
  TraceResult& operator=(TraceResult&&) = default;

  // Allow copies explicitly using the constructor.
  explicit TraceResult(const TraceResult&) = default;
  TraceResult& operator=(const TraceResult&) = delete;

  constexpr auto operator<=>(const TraceResult&) const = default;

  Secondary s;
  Ctds ctd;  // May be empty if CTDs were not computed.
};

}  // namespace mrna::trace

#endif  // API_TRACE_TRACE_H_
