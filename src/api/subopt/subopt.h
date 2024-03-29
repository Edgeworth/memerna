// Copyright 2016 E.
#ifndef API_SUBOPT_SUBOPT_H_
#define API_SUBOPT_SUBOPT_H_

#include <functional>
#include <utility>

#include "api/trace/trace.h"
#include "model/energy.h"

namespace mrna::subopt {

struct SuboptResult {
  SuboptResult() = default;
  ~SuboptResult() = default;
  SuboptResult(Energy energy, trace::TraceResult tb) : energy(energy), tb(std::move(tb)) {}

  SuboptResult(SuboptResult&&) = default;
  SuboptResult& operator=(SuboptResult&&) = default;

  // Allow copies explicitly using the constructor.
  explicit SuboptResult(const SuboptResult&) = default;
  SuboptResult& operator=(const SuboptResult&) = delete;

  constexpr auto operator<=>(const SuboptResult&) const = default;

  Energy energy{};  // Put this first so naive sort is by energy.
  trace::TraceResult tb;
};

using SuboptCallback = std::function<void(const SuboptResult&)>;

}  // namespace mrna::subopt

#endif  // API_SUBOPT_SUBOPT_H_
