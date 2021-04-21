// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#ifndef MEMERNA_FOLD_CONSTANTS_H
#define MEMERNA_FOLD_CONSTANTS_H

namespace memerna {
namespace fold {
namespace internal {

// DP arrays
enum : int8_t {
  DP_P,     // For the paired array.
  DP_U,     // For the unpaired array.
  DP_U2,    // Contains at least two branches.
  DP_U_WC,  // Unpaired but must start with a branch not in a CTD that is not GU.
  DP_U_GU,  // Unpaired but must start with a branch not in a CTD that is GU.
  DP_U_RCOAX,  // Unpaired but must start with a branch in a RCOAX - includes energy for it.
  DP_SIZE
};

enum : int8_t {
  EXT,
  EXT_WC,     // Must start with a branch not involved in an interaction that is Watson-Crick
  EXT_GU,     // Must start with a branch not involved in an interaction that is GU
  EXT_RCOAX,  // Must start with a branch, that branch is involved backwards in a RCOAX stack.
  EXT_SIZE
};

// Split candidates up into several lists.
// In general, for each array we need a new candidate list (except for U and U2 which mirror each
// other very closely). We also need another candidate list for forward RCOAX since we can't use
// its energy value directly, not knowing it. Same with flush coaxial stacks.
enum : int8_t {
  CAND_P_MISMATCH,  // No monotonicity.
  CAND_P_OUTER,     // No monotonicity.
  CAND_P_FLUSH,     // No monotonicity.
  CAND_U,
  CAND_U_LCOAX,      // No monotonicity.
  CAND_U_RCOAX_FWD,  // No monotonicity.
  CAND_U_WC_FLUSH,   // No monotonicity.
  CAND_U_GU_FLUSH,   // No monotonicity.
  CAND_U_WC,
  CAND_U_GU,
  CAND_U_RCOAX,
  CAND_SIZE
};

enum : int8_t {
  CAND_EN_P_MISMATCH, CAND_EN_P_OUTER, CAND_EN_P_FLUSH, CAND_EN_SIZE
};

const int MAX_STRUCTURES = std::numeric_limits<int>::max() / 4;

}
}
}

#endif  // MEMERNA_FOLD_CONSTANTS_H
