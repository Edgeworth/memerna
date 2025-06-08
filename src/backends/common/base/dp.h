// Copyright 2025 Eliot Courtney.
#ifndef BACKENDS_COMMON_BASE_DP_H_
#define BACKENDS_COMMON_BASE_DP_H_

#include <algorithm>

#include "api/trace/trace.h"
#include "util/array.h"

namespace mrna::md::base {

using trace::TraceResult;

// DP arrays
enum : int8_t {
  DP_P,  // For the paired array.
  DP_U,  // For the unpaired array. Contains at least one branch. First and last base may be paired.
  DP_U2,  // Contains at least two branches.
  DP_U_WC,  // Unpaired but must start with a branch not in a CTD that is not GU.
  DP_U_GU,  // Unpaired but must start with a branch not in a CTD that is GU.
  DP_U_RC,  // Unpaired but must start with a branch in a right coaxial stack - includes energy
            // for it.
  DP_SIZE
};

using DpArray = Array2D1S<Energy, DP_SIZE>;

enum : int8_t {
  EXT,
  EXT_WC,  // Must start with a branch not involved in an interaction that is Watson-Crick
  EXT_GU,  // Must start with a branch not involved in an interaction that is GU
  EXT_RC,  // Must start with a branch involved backwards in a right coaxial stack.
  EXT_SIZE
};

using ExtArray = Array1D1S<Energy, EXT_SIZE>;

// Split candidates up into several lists. In general, for each array we need a
// new candidate list (except for U and U2 which mirror each other very
// closely). We also need another candidate list for forward right coaxial since
// we can't use its energy value directly, not knowing it. Same with flush
// coaxial stacks. |...| denotes which area is replaceable.
enum : int8_t {
  // Monotonic, because we place Unpaired on the right.
  // Replaceability: Unpaired
  // Energy: Fully determined, since we have all the bases available.
  // Left inner coax: Mismatch mediated coaxial stack on the left with the outer loop.
  // (|.(   ).|   )
  CAND_P_LIC,

  // Monotonic, because we place Unpaired on the right.
  // Replaceability: Unpaired
  // Energy: Partial, since we don't know the bases on the right.
  // Left outer coax: Outer loop outer mismatch mediated coaxial stack with a branch on the left
  // (|.(   )|   .)
  CAND_P_LOC,

  // Monotonic, because we place Unpaired on the right.
  // Replaceability: Unpaired
  // Energy: Partial, since we don't know the bases on the right.
  // Left flush coax: Outer loop forms a flush coaxial stack with a branch on the left
  // (|(   )|   )
  CAND_P_LFC,

  // Monotonic, because we place Unpaired on the right.
  // Replaceability: NoCTD, Dangle3', Dangle5', Mismatch
  // Energy: Fully determined, since just variations on Paired.
  // Unpaired
  // |<   >|<   >
  CAND_U,

  // Not monotonic, because we are forced to place a branch on the right.
  // Replaceability: Unpaired
  // Energy: Fully determined, since we have all the bases available.
  // Left coax
  // |.(   ).|<(   ) >
  CAND_U_LC,

  // Not monotonic, because we are forced to place a branch on the right.
  // Replaceability: Unpaired
  // Energy: Partial, since we don't know the bases on the right.
  // Right coax forward
  // |(   )|<.(   ). >
  CAND_U_RC_FWD,

  // Not monotonic, because we are forced to place a branch on the right.
  // Replaceability: Unpaired
  // Energy: Fully determined, since we know it is a WC pair on the right.
  // Flush coax with WC branch
  // |(   )|(<   ) >
  CAND_U_LFC_WC,

  // Not monotonic, because we are forced to place a branch on the right.
  // Replaceability: Unpaired
  // Energy: Fully determined, since we know it is a GU pair on the right.
  // Flush coax with GU branch
  // |(   )|(<   ) >
  CAND_U_LFC_GU,

  // Monotonic, because we place Unpaired on the right.
  // Replaceability: Only itself.
  // Energy: Fully determined, since just variations on Paired with no CTDs.
  // UnpairedWC
  // |<   >|<   >
  CAND_U_WC,

  // Monotonic, because we place Unpaired on the right.
  // Replaceability: Only itself.
  // Energy: Fully determined, since just variations on Paired with no CTDs.
  // UnpairedGU
  // |(   )|<   >
  CAND_U_GU,

  // Monotonic, because we place Unpaired on the right.
  // Replaceability: Only itself.
  // Energy: Fully determined, because we know all the bases.
  // UnpairedRcoax
  // (   )|<.( * ). >|
  CAND_U_RC,

  CAND_SIZE
};

enum : int8_t {
  // Monotonic, because we place Unpaired on the left.
  // Replaceability: Unpaired
  // Energy: Fully determined, since we have all the bases available.
  // Right inner coax: Mismatch mediated coaxial stack on the right with the outer loop.
  // (   |.(   ).|)
  CAND_EN_P_RIC,

  // Monotonic, because we place Unpaired on the left.
  // Replaceability: Unpaired
  // Energy: Partial, since we don't know the bases on the left.
  // Right outer coax: Outer loop outer mismatch mediated coaxial stack with a branch on the right
  // (.   |(   ).|)
  CAND_EN_P_ROC,

  // Monotonic, because we place Unpaired on the left.
  // Replaceability: Unpaired
  // Energy: Partial, since we don't know the bases on the left.
  // Right flush coax: Outer loop forms a flush coaxial stack with a branch on the right
  // (   |(   )|)
  CAND_EN_P_RFC,

  CAND_EN_SIZE
};

// Index into the DP tables. `en` set to -1 is used to indicate using the external table. `st` can
// be N on the external loop.
struct DpIndex {
  Index st{-1}, en{-1}, a{-1};

  DpIndex() = default;
  DpIndex(int st_, int en_, int a_) : st(Index(st_)), en(Index(en_)), a(Index(a_)) {
    assert(st_ == st && en_ == en && a == a_);
  }

  constexpr auto operator<=>(const DpIndex&) const = default;

  [[nodiscard]] constexpr std::size_t LinearIndex(std::size_t n) const {
    assert(a >= 0 && a < int(MaxArrayCount()));
    assert(st >= 0 && st <= int(n));
    assert(en >= -1 && en < int(n));
    return a + MaxArrayCount() * (en + 1) + MaxArrayCount() * (n + 1) * st;
  }

  [[nodiscard]] constexpr static std::size_t MaxLinearIndex(std::size_t n) {
    return MaxArrayCount() + MaxArrayCount() * (n + 1) + MaxArrayCount() * (n + 1) * (n + 1);
  }

 private:
  [[nodiscard]] constexpr static std::size_t MaxArrayCount() {
    return std::max(static_cast<std::size_t>(DP_SIZE), static_cast<std::size_t>(EXT_SIZE));
  }
};

struct Cand {
  Energy energy;
  int idx;
};

struct DpState {
  DpArray dp;
  ExtArray ext;

  [[nodiscard]] constexpr Energy Index(const DpIndex& idx) const {
    if (idx.en == -1) {
      return ext[idx.st][idx.a];
    }
    return dp[idx.st][idx.en][idx.a];
  }
};

struct Expansion {
  // Extra energy of this expansion compared to the best choice.
  Energy delta = {ZERO_E};

  // st == -1 used to mean none - using optional here is like a 40% perf hit.
  DpIndex idx0{};
  DpIndex idx1{};
  IndexCtd ctd0{};
  IndexCtd ctd1{};

  bool operator<(const Expansion& o) const { return delta < o.delta; }
};

// DP arrays
enum : int8_t { PT_P, PT_U, PT_U2, PT_U_WC, PT_U_GU, PT_U_RC, PT_SIZE };

using BoltzDpArray = Array2D1S<BoltzEnergy, PT_SIZE>;

enum : int8_t {
  PTEXT_R,
  PTEXT_L,
  PTEXT_R_WC,  // Must start with a branch not involved in an interaction that is Watson-Crick
  PTEXT_R_GU,  // Must start with a branch not involved in an interaction that is GU
  // Must start with a branch, that branch is involved backwards in a right coaxial stack.
  PTEXT_R_RC,
  PTEXT_L_WC,
  PTEXT_L_GU,
  PTEXT_L_LCOAX,
  PTEXT_SIZE
};

using BoltzExtArray = Array1D1S<BoltzEnergy, PTEXT_SIZE>;

struct PfnState {
  BoltzDpArray dp;
  BoltzExtArray ext;
};

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_DP_H_
