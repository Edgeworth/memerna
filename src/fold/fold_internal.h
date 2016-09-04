#ifndef MEMERNA_FOLD_INTERNAL_H
#define MEMERNA_FOLD_INTERNAL_H

#include "common.h"
#include "constants.h"
#include "array.h"
#include "energy/energy.h"
#include "energy/energy_model.h"

namespace memerna {
namespace fold {
namespace internal {

inline bool ViableFoldingPair(const primary_t& r, int st, int en) {
  return CanPair(r[st], r[en]) &&
      ((en - st - 3 >= constants::HAIRPIN_MIN_SZ && CanPair(r[st + 1], r[en - 1])) ||
          (st > 0 && en < int(r.size() - 1) && CanPair(r[st - 1], r[en + 1])));
}

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

struct hairpin_precomp_t {
  static const int MAX_SPECIAL_HAIRPIN_SZ = 6;
  hairpin_precomp_t() : num_c(0) {
    memset(special, constants::MAX_E & 0xFF, sizeof(special));
  }

  energy_t special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c;
};

struct precomp_t {
  energy_t augubranch[4][4];
  energy_t min_mismatch_coax;
  energy_t min_flush_coax;
  energy_t min_twoloop_not_stack;

  std::vector<hairpin_precomp_t> hairpin;
};

int MaxNumContiguous(const primary_t& r);
precomp_t PrecomputeData(const primary_t& r, const energy::EnergyModel& em);

}
}
}

#endif //MEMERNA_FOLD_INTERNAL_H
