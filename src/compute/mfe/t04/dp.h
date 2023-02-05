// Copyright 2021 Eliot Courtney.
#ifndef COMPUTE_MFE_T04_DP_H_
#define COMPUTE_MFE_T04_DP_H_

#include <cassert>

#include "model/base.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "util/array.h"

namespace mrna::mfe::t04 {

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

// Index into the DP tables.
// Use int16_t here to save memory.
struct Index {
  int16_t st{-1}, en{-1}, a{-1};

  Index() = default;
  Index(int st_, int en_, int a_) : st(int16_t(st_)), en(int16_t(en_)), a(int16_t(a_)) {
    assert(st_ == st && en_ == en && a == a_);
  }

  constexpr auto operator<=>(const Index&) const = default;
};

// Describes a CTD at a particular index.
struct IndexCtd {
  IndexCtd() = default;
  IndexCtd(int idx_, Ctd ctd_) : idx(int16_t(idx_)), ctd(ctd_) { assert(idx_ == idx); }

  int16_t idx{-1};
  Ctd ctd{CTD_NA};
};

struct DpState {
  DpArray dp;
  ExtArray ext;
};

}  // namespace mrna::mfe::t04

#endif  // COMPUTE_MFE_T04_DP_H_
