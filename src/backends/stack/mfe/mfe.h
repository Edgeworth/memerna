// Copyright 2023 Eliot Courtney.

#ifndef BACKENDS_STACK_MFE_MFE_H_
#define BACKENDS_STACK_MFE_MFE_H_

#include <variant>

#include "backends/common/base/dp.h"
#include "backends/stack/energy/model.h"
#include "model/primary.h"
#include "util/util.h"

namespace mrna::md::stack {

struct DpState {
  // base state is reused.
  base::DpState base;
  // DP for holding the best answer given (st, en) is paired but we can't
  // start a stack of length 2 or more from here.
  Array2D<Energy> nostack;
  // DP for computing the best stack of a given length.
  Array3D<Energy> penult;
};

// Use int16_t here to save memory.
struct PenultimateIndex {
  Index st, en, len;

  PenultimateIndex(int st_, int en_, int len_) : st(Index(st_)), en(Index(en_)), len(Index(len_)) {
    assert(st_ == st && en_ == en && len == len_);
  }

  constexpr auto operator<=>(const PenultimateIndex&) const = default;

  [[nodiscard]] constexpr std::size_t LinearIndex(std::size_t n) const {
    assert(st >= 0 && st <= int(n));
    assert(en >= 0 && en <= int(n));
    assert(len >= 0 && len <= int(n));

    return len + (n + 1) * en + (n + 1) * (n + 1) * st;
  }

  [[nodiscard]] constexpr static std::size_t MaxLinearIndex(std::size_t n) {
    return (n + 1) * (n + 1) * (n + 1);
  }
};

struct NoStackIndex {
  Index st, en;

  NoStackIndex(int st_, int en_) : st(Index(st_)), en(Index(en_)) {
    assert(st_ == st && en_ == en);
  }

  constexpr auto operator<=>(const NoStackIndex&) const = default;

  [[nodiscard]] constexpr std::size_t LinearIndex(std::size_t n) const {
    assert(st >= 0 && st <= int(n));
    assert(en >= 0 && en <= int(n));
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

void MfeDebug(const Primary& r, const Model::Ptr& m, DpState& state);

Energy MfeExterior(const Primary& r, const Model::Ptr& m, DpState& state);

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_MFE_MFE_H_
