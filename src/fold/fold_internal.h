#ifndef MEMERNA_FOLD_INTERNAL_H
#define MEMERNA_FOLD_INTERNAL_H

#include "common.h"
#include "constants.h"
#include "array.h"
#include "energy/energy.h"
#include "fold/fold.h"

namespace memerna {
namespace fold {
namespace internal {

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
  CAND_U_WC_FLUSH,  // No monotonicity.
  CAND_U_GU_FLUSH,  // No monotonicity.
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

const int MAX_SPECIAL_HAIRPIN_SZ = 6;

struct hairpin_precomp_t {
  hairpin_precomp_t() : num_c(0) {
    memset(special, constants::MAX_E & 0xFF, sizeof(special));
  }

  energy_t special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c;
};

array3d_t<energy_t, DP_SIZE> ComputeTables3(const primary_t& r);
array3d_t<energy_t, DP_SIZE> ComputeTables2(const primary_t& r);
array3d_t<energy_t, DP_SIZE> ComputeTables1(const primary_t& r);
array3d_t<energy_t, DP_SIZE> ComputeTables0(const primary_t& r);

int MaxNumContiguous(const primary_t& r);
void InitFold(const primary_t& r);

energy_t FastTwoLoop(const primary_t& r, int ost, int oen, int ist, int ien);

std::vector<hairpin_precomp_t> PrecomputeFastHairpin(const primary_t& r);
energy_t FastHairpin(const primary_t& r, int st, int en, const std::vector<hairpin_precomp_t>& precomp);

inline bool ViableFoldingPair(const primary_t& r, int st, int en) {
  return CanPair(r[st], r[en]) &&
      ((en - st - 3 >= constants::HAIRPIN_MIN_SZ && CanPair(r[st + 1], r[en - 1])) ||
          (st > 0 && en < int(r.size() - 1) && CanPair(r[st - 1], r[en + 1])));
}

}
}
}

#endif //MEMERNA_FOLD_INTERNAL_H
