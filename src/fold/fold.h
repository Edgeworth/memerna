#ifndef MEMERNA_FOLD_H
#define MEMERNA_FOLD_H

#include "constants.h"
#include "common.h"
#include "energy/energy.h"

namespace memerna {
namespace fold {

const int MAX_SPECIAL_HAIRPIN_SZ = 6;

// Non continuous -2.1, -4 for WC, -16 for terminal mismatch.
// TODO: compute in InitFold
const energy_t MIN_MISMATCH_COAX = -21 - 4 - 16;
const energy_t MIN_FLUSH_COAX = -34;

// DP arrays
enum {
  DP_P,  // For the paired array.
  DP_U,  // For the unpaired array.
  DP_U2, // Contains at least two branches.
  DP_U_WC,  // Unpaired but must start with a branch not involved in a CTD interaction that is not GU.
  DP_U_GU,  // Unpaired but must start with a branch not involved in a CTD interaction that is GU.
  DP_U_RCOAX,  // Unpaired but must start with a branch involved in a right coaxial stack - includes energy for it.
  DP_SIZE
};

enum {
  EXT,
  EXT_WC,
  EXT_GU,
  EXT_RCOAX,
  EXT_SIZE
};

// Split candidates up into several lists.
// In general, for each array we need a new candidate list (except for U and U2 which mirror each other very
// closely). We also need another candidate list for forward RCOAX since we can't use its energy value directly, not
// knowing it. Same with flush coaxial stacks.
enum {
  CAND_P_MISMATCH,  // No monotonicity.
  CAND_P_OUTER,  // No monotonicity.
  CAND_P_FLUSH,  // No monotonicity.
  CAND_U,
  CAND_U_LCOAX,  // No monotonicity.
  CAND_U_RCOAX_FWD,  // No monotonicity.
  CAND_U_FLUSH,  // No monotonicity.
  CAND_U_WC,
  CAND_U_GU,
  CAND_U_RCOAX,
  CAND_SIZE
};

enum {
  CAND_EN_P_MISMATCH,
  CAND_EN_P_OUTER,
  CAND_EN_P_FLUSH,
  CAND_EN_SIZE
};

struct cand_t {
  energy_t energy;
  int idx;
};


void InitFold();
energy_t FastTwoLoop(int ost, int oen, int ist, int ien);

struct hairpin_precomp_t {
  hairpin_precomp_t() : num_c(0) {
    memset(special, constants::MAX_E & 0xFF, sizeof(special));
  }
  energy_t special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c;
};

std::vector<hairpin_precomp_t> PrecomputeFastHairpin();
energy_t FastHairpin(int st, int en, const std::vector<hairpin_precomp_t>& precomp);


energy_t Fold();
energy_t Fold3();
energy_t Fold2();
energy_t Fold1();
energy_t FoldSlow();
energy_t FoldBruteForce();

inline energy_t Fold(const rna_t& rna) {
  SetRna(rna);
  return Fold();
}

inline bool IsNotLonely(int st, int en) {
  return (en - st - 3 >= constants::HAIRPIN_MIN_SZ && CanPair(r[st + 1], r[en - 1])) ||
      (st > 0 && en < int(r.size() - 1) && CanPair(r[st - 1], r[en + 1]));
}

}
}

#endif //MEMERNA_FOLD_H
