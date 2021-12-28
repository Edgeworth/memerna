// Copyright 2021 E.
#ifndef COMPUTE_DP_H_
#define COMPUTE_DP_H_

#include <cassert>

#include "model/base.h"
#include "model/globals.h"

namespace mrna {

// DP arrays
enum : int8_t {
  DP_P,  // For the paired array.
  DP_U,  // For the unpaired array.
  DP_U2,  // Contains at least two branches.
  DP_U_WC,  // Unpaired but must start with a branch not in a CTD that is not GU.
  DP_U_GU,  // Unpaired but must start with a branch not in a CTD that is GU.
  DP_U_RCOAX,  // Unpaired but must start with a branch in a RCOAX - includes energy for it.
  DP_SIZE
};

enum : int8_t {
  EXT,
  EXT_WC,  // Must start with a branch not involved in an interaction that is Watson-Crick
  EXT_GU,  // Must start with a branch not involved in an interaction that is GU
  EXT_RCOAX,  // Must start with a branch, that branch is involved backwards in a RCOAX stack.
  EXT_SIZE
};

// Split candidates up into several lists.
// In general, for each array we need a new candidate list (except for U and U2 which mirror each
// other very closely). We also need another candidate list for forward RCOAX since we can't use
// its energy value directly, not knowing it. Same with flush coaxial stacks.
enum : int8_t {
  CAND_P_MISMATCH,  // No monotonicity.
  CAND_P_OUTER,  // No monotonicity.
  CAND_P_FLUSH,  // No monotonicity.
  CAND_U,
  CAND_U_LCOAX,  // No monotonicity.
  CAND_U_RCOAX_FWD,  // No monotonicity.
  CAND_U_WC_FLUSH,  // No monotonicity.
  CAND_U_GU_FLUSH,  // No monotonicity.
  CAND_U_WC,
  CAND_U_GU,
  CAND_U_RCOAX,
  CAND_SIZE
};

enum : int8_t { CAND_EN_P_MISMATCH, CAND_EN_P_OUTER, CAND_EN_P_FLUSH, CAND_EN_SIZE };

// Index into the DP tables.
// Use int16_t here to save memory.
struct index_t {
  int16_t st, en, a;

  index_t() : st(-1), en(-1), a(-1) {}
  index_t(int st_, int en_, int a_) : st(int16_t(st_)), en(int16_t(en_)), a(int16_t(a_)) {
    assert(st_ == st && en_ == en && a == a_);
  }

  bool operator==(const index_t& o) const { return st == o.st && en == o.en && a == o.a; }
  bool operator!=(const index_t& o) const { return !(*this == o); }
  bool operator<(const index_t& o) const {
    if (st != o.st) return st < o.st;
    if (en != o.en) return en < o.en;
    return a < o.a;
  }
};

// Describes a CTD at a particular index.
struct ctd_idx_t {
  ctd_idx_t() : idx(-1), ctd(CTD_NA) {}
  ctd_idx_t(int idx_, Ctd ctd_) : idx(int16_t(idx_)), ctd(ctd_) { assert(idx_ == idx); }

  int16_t idx;
  Ctd ctd;
};

}  // namespace mrna

#endif  // COMPUTE_DP_H_