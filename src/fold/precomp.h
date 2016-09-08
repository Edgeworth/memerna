#ifndef MEMERNA_PRECOMP_H
#define MEMERNA_PRECOMP_H

#include "common.h"
#include "constants.h"
#include "energy/energy_model.h"

namespace memerna {
namespace fold {
namespace internal {

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
// Must have global state set.
energy_t FastTwoLoop(int ost, int oen, int ist, int ien);
energy_t FastHairpin(int st, int en);

}
}
}

#endif  // MEMERNA_PRECOMP_H
