// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_PARTITION_PARTITION_H_
#define COMPUTE_PARTITION_PARTITION_H_

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "compute/energy/fast_energy.h"
#include "compute/energy/model.h"
#include "util/array.h"

namespace mrna::partition {

enum : int8_t { PT_P, PT_U, PT_U2, PT_U_WC, PT_U_GU, PT_U_RCOAX, PT_SIZE };

enum : int8_t {
  PTEXT_R,
  PTEXT_L,
  PTEXT_R_WC,  // Must start with a branch not involved in an interaction that is Watson-Crick
  PTEXT_R_GU,  // Must start with a branch not involved in an interaction that is GU
  PTEXT_R_RCOAX,  // Must start with a branch, that branch is involved backwards in a RCOAX stack.
  PTEXT_L_WC,
  PTEXT_L_GU,
  PTEXT_L_LCOAX,
  PTEXT_SIZE
};

typedef Array3D<PEnergy, 1> Probabilities;

struct Partition {
  Array3D<PEnergy, 1> p;
  PEnergy q;
};

Probabilities ComputeProbabilities(const Partition& partition);

namespace internal {

struct PEnergyModel {
  PEnergy stack[4][4][4][4];
  PEnergy terminal[4][4][4][4];
  PEnergy internal_init[energy::EnergyModel::INITIATION_CACHE_SZ];
  PEnergy internal_1x1[4][4][4][4][4][4];
  PEnergy internal_1x2[4][4][4][4][4][4][4];
  PEnergy internal_2x2[4][4][4][4][4][4][4][4];
  PEnergy internal_2x3_mismatch[4][4][4][4];
  PEnergy internal_other_mismatch[4][4][4][4];
  PEnergy internal_asym;
  PEnergy internal_augu_penalty;
  PEnergy bulge_init[energy::EnergyModel::INITIATION_CACHE_SZ];
  PEnergy bulge_special_c;
  PEnergy hairpin_init[energy::EnergyModel::INITIATION_CACHE_SZ];
  PEnergy hairpin_uu_ga_first_mismatch, hairpin_gg_first_mismatch, hairpin_special_gu_closure,
      hairpin_c3_loop, hairpin_all_c_a, hairpin_all_c_b;
  std::unordered_map<std::string, PEnergy> hairpin;
  PEnergy multiloop_hack_a, multiloop_hack_b;
  PEnergy dangle5[4][4][4];
  PEnergy dangle3[4][4][4];
  PEnergy coax_mismatch_non_contiguous, coax_mismatch_wc_bonus, coax_mismatch_gu_bonus;
  PEnergy augu_penalty;

  PEnergyModel() = default;
  explicit PEnergyModel(const energy::EnergyModel& em);

  PEnergy InternalLoopAuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    return IsAuGu(stb, enb) ? internal_augu_penalty : 1.0;
  }

  PEnergy AuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    return IsAuGu(stb, enb) ? augu_penalty : 1.0;
  }

  PEnergy MismatchCoaxial(
      Base five_top, Base mismatch_top, Base mismatch_bot, Base three_bot) const {
    assert(IsBase(five_top) && IsBase(mismatch_top) && IsBase(mismatch_bot) && IsBase(three_bot));
    PEnergy coax =
        terminal[five_top][mismatch_top][mismatch_bot][three_bot] * coax_mismatch_non_contiguous;
    if (IsWatsonCrick(mismatch_top, mismatch_bot))
      coax *= coax_mismatch_wc_bonus;
    else if (IsGu(mismatch_top, mismatch_bot))
      coax *= coax_mismatch_gu_bonus;
    return coax;
  }
};

struct Precomp {
  PEnergy augubranch[4][4];
  PEnergyModel em;
  std::vector<energy::HairpinPrecomp<PEnergy, -1>> hairpin;
};

Precomp PrecomputeData(const Primary& r, const energy::EnergyModel& em);
PEnergy FastHairpin(int st, int en);
PEnergy FastTwoLoop(int ost, int oen, int ist, int ien);
void Partition1();
void Partition0();
void Exterior();

// Only works with [0, 2N).
inline int FastMod(int a, int m) {
  if (a >= m) return a - m;
  return a;
}

}  // namespace internal

}  // namespace mrna::partition

#endif  // COMPUTE_PARTITION_PARTITION_H_
