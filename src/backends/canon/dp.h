// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_DP_H_
#define BACKENDS_CANON_DP_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "model/base.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "util/array.h"
#include "util/util.h"

namespace mrna::md::canon {

using DpArrayId = uint8_t;
constexpr DpArrayId INVALID_DP_ARRAY_ID = std::numeric_limits<DpArrayId>::max();

enum : DpArrayId {
  DP_P,
  DP_U,
  DP_U2,
  DP_CX,
  DP_SIZE,
};

using DpArray = Array2D1S<Energy, DP_SIZE>;

enum : DpArrayId {
  EXT,
  EXT_SIZE,
};

using ExtArray = Array1D1S<Energy, EXT_SIZE>;

struct DpIndex {
  Index st{INVALID_INDEX}, en{INVALID_INDEX};
  DpArrayId a{INVALID_DP_ARRAY_ID};

  constexpr DpIndex() = default;
  constexpr DpIndex(int st_, int en_, DpArrayId a_)
      : st(As<Index>(st_)), en(As<Index>(en_)), a(a_) {}

  constexpr auto operator<=>(const DpIndex&) const = default;

  [[nodiscard]] constexpr bool IsValid() const { return st != INVALID_INDEX; }

  [[nodiscard]] constexpr bool IsExternal() const { return en == INVALID_INDEX; }

  [[nodiscard]] constexpr std::size_t LinearIndex(std::size_t n) const {
    assert(a < MaxArrayCount());
    assert(st >= 0 && std::size_t(st) <= n);
    assert(IsExternal() || (en >= 0 && std::size_t(en) < n));
    const std::size_t en_idx = IsExternal() ? 0 : en + 1;
    return a + MaxArrayCount() * en_idx + MaxArrayCount() * (n + 1) * st;
  }

  [[nodiscard]] constexpr static std::size_t MaxLinearIndex(std::size_t n) {
    return MaxArrayCount() * (n + 1) * (n + 1);
  }

 private:
  [[nodiscard]] constexpr static std::size_t MaxArrayCount() {
    return std::max(static_cast<std::size_t>(DP_SIZE), static_cast<std::size_t>(EXT_SIZE));
  }
};

struct DpState {
  DpArray dp;
  ExtArray ext;

  [[nodiscard]] constexpr Energy Index(const DpIndex& idx) const {
    if (idx.IsExternal()) return ext[idx.st][idx.a];
    return dp[idx.st][idx.en][idx.a];
  }
};

struct Expansion {
  Energy delta = ZERO_E;
  DpIndex idx0{};
  DpIndex idx1{};
  DpIndex idx2{};
  IndexCtd ctd0{};
  IndexCtd ctd1{};
  IndexCtd ctd2{};

  bool operator<(const Expansion& o) const { return delta < o.delta; }
};

enum : DpArrayId {
  PT_P,
  PT_U,
  PT_U2,
  PT_CX,
  PT_SIZE,
};

using BoltzDpArray = Array2D1S<BoltzEnergy, PT_SIZE>;

enum : DpArrayId {
  PTEXT_R,
  PTEXT_L,
  PTEXT_SIZE,
};

using BoltzExtArray = Array1D1S<BoltzEnergy, PTEXT_SIZE>;

struct PfnState {
  BoltzDpArray dp;
  BoltzExtArray ext;
};

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_DP_H_
