// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_SUBOPT_SUBOPT_H_
#define COMPUTE_SUBOPT_SUBOPT_H_

#include <functional>
#include <utility>

#include "compute/traceback/traceback.h"
#include "model/energy.h"

namespace mrna::subopt {

struct SuboptResult {
  SuboptResult() = default;
  ~SuboptResult() = default;
  SuboptResult(Energy energy, tb::TraceResult tb) : energy(energy), tb(std::move(tb)) {}

  SuboptResult(SuboptResult&&) = default;
  SuboptResult& operator=(SuboptResult&&) = default;

  // Allow copies explicitly using the constructor.
  explicit SuboptResult(const SuboptResult&) = default;
  SuboptResult& operator=(const SuboptResult&) = delete;

  constexpr auto operator<=>(const SuboptResult&) const = default;

  Energy energy{};  // Put this first so naive sort is by energy.
  tb::TraceResult tb;
};

using SuboptCallback = std::function<void(const SuboptResult&)>;

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_SUBOPT_H_
