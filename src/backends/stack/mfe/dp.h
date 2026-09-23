// Copyright 2026 Eliot Courtney.

#ifndef BACKENDS_STACK_MFE_DP_H_
#define BACKENDS_STACK_MFE_DP_H_

#include <cassert>
#include <cstddef>
#include <variant>

#include "backends/common/base/dp.h"
#include "model/base.h"
#include "model/energy.h"
#include "util/array.h"
#include "util/util.h"

namespace mrna::md::stack {

struct DpState {
  // base state is reused.
  base::DpState base;
  // DP for holding the best answer given (st, en) is paired but we can't
  // start a stack of length 2 or more from here.
  Array2D<Energy> nostack;
  // DP for the best stack of length 2 or more starting at (st, en). Includes the
  // penultimate penalty at the inner end of the stack, but not the one at
  // (st, en), which only applies if the stack does not continue outwards.
  Array2D<Energy> penult;
};

struct PenultimateIndex {
  Index st, en;

  PenultimateIndex(int st_, int en_) : st(As<Index>(st_)), en(As<Index>(en_)) {}

  constexpr auto operator<=>(const PenultimateIndex&) const = default;

  [[nodiscard]] constexpr std::size_t LinearIndex(std::size_t n) const {
    assert(st >= 0 && std::size_t(st) <= n);
    assert(en >= 0 && std::size_t(en) <= n);
    return en + (n + 1) * st;
  }

  [[nodiscard]] constexpr static std::size_t MaxLinearIndex(std::size_t n) {
    return (n + 1) * (n + 1);
  }
};

struct NoStackIndex {
  Index st, en;

  NoStackIndex(int st_, int en_) : st(As<Index>(st_)), en(As<Index>(en_)) {}

  constexpr auto operator<=>(const NoStackIndex&) const = default;

  [[nodiscard]] constexpr std::size_t LinearIndex(std::size_t n) const {
    assert(st >= 0 && std::size_t(st) <= n);
    assert(en >= 0 && std::size_t(en) <= n);
    return en + (n + 1) * st;
  }

  [[nodiscard]] constexpr static std::size_t MaxLinearIndex(std::size_t n) {
    return (n + 1) * (n + 1);
  }
};

// Index into the DP tables.
using DpIndex = std::variant<base::DpIndex, NoStackIndex, PenultimateIndex>;

constexpr std::size_t LinearIndex(const DpIndex& index, std::size_t n) {
  return std::visit(overloaded{
                        [n](const base::DpIndex& idx) { return idx.LinearIndex(n); },
                        [n](const NoStackIndex& idx) {
                          return idx.LinearIndex(n) + base::DpIndex::MaxLinearIndex(n);
                        },
                        [n](const PenultimateIndex& idx) {
                          return idx.LinearIndex(n) + base::DpIndex::MaxLinearIndex(n) +
                              NoStackIndex::MaxLinearIndex(n);
                        },
                    },
      index);
}

constexpr std::size_t MaxLinearIndex(std::size_t n) {
  return base::DpIndex::MaxLinearIndex(n) + NoStackIndex::MaxLinearIndex(n) +
      PenultimateIndex::MaxLinearIndex(n);
}

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_MFE_DP_H_
