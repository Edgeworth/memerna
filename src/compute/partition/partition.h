// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_PARTITION_PARTITION_H_
#define COMPUTE_PARTITION_PARTITION_H_

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "common.h"
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

typedef array3d_t<penergy_t, 1> probabilities_t;

struct partition_t {
  array3d_t<penergy_t, 1> p;
  penergy_t q;
};

probabilities_t ComputeProbabilities(const partition_t& partition);

namespace internal {

struct penergy_model_t {
  penergy_t stack[4][4][4][4];
  penergy_t terminal[4][4][4][4];
  penergy_t internal_init[energy::EnergyModel::INITIATION_CACHE_SZ];
  penergy_t internal_1x1[4][4][4][4][4][4];
  penergy_t internal_1x2[4][4][4][4][4][4][4];
  penergy_t internal_2x2[4][4][4][4][4][4][4][4];
  penergy_t internal_2x3_mismatch[4][4][4][4];
  penergy_t internal_other_mismatch[4][4][4][4];
  penergy_t internal_asym;
  penergy_t internal_augu_penalty;
  penergy_t bulge_init[energy::EnergyModel::INITIATION_CACHE_SZ];
  penergy_t bulge_special_c;
  penergy_t hairpin_init[energy::EnergyModel::INITIATION_CACHE_SZ];
  penergy_t hairpin_uu_ga_first_mismatch, hairpin_gg_first_mismatch, hairpin_special_gu_closure,
      hairpin_c3_loop, hairpin_all_c_a, hairpin_all_c_b;
  std::unordered_map<std::string, penergy_t> hairpin;
  penergy_t multiloop_hack_a, multiloop_hack_b;
  penergy_t dangle5[4][4][4];
  penergy_t dangle3[4][4][4];
  penergy_t coax_mismatch_non_contiguous, coax_mismatch_wc_bonus, coax_mismatch_gu_bonus;
  penergy_t augu_penalty;

  penergy_model_t() = default;
  explicit penergy_model_t(const energy::EnergyModel& em);

  penergy_t InternalLoopAuGuPenalty(base_t stb, base_t enb) const {
    assert(IsBase(stb) && IsBase(enb));
    return IsAuGu(stb, enb) ? internal_augu_penalty : 1.0;
  }

  penergy_t AuGuPenalty(base_t stb, base_t enb) const {
    assert(IsBase(stb) && IsBase(enb));
    return IsAuGu(stb, enb) ? augu_penalty : 1.0;
  }

  penergy_t MismatchCoaxial(
      base_t five_top, base_t mismatch_top, base_t mismatch_bot, base_t three_bot) const {
    assert(IsBase(five_top) && IsBase(mismatch_top) && IsBase(mismatch_bot) && IsBase(three_bot));
    penergy_t coax =
        terminal[five_top][mismatch_top][mismatch_bot][three_bot] * coax_mismatch_non_contiguous;
    if (IsWatsonCrick(mismatch_top, mismatch_bot))
      coax *= coax_mismatch_wc_bonus;
    else if (IsGu(mismatch_top, mismatch_bot))
      coax *= coax_mismatch_gu_bonus;
    return coax;
  }
};

struct precomp_t {
  penergy_t augubranch[4][4];
  penergy_model_t em;
  std::vector<energy::hairpin_precomp_t<penergy_t, -1>> hairpin;
};

precomp_t PrecomputeData(const primary_t& r, const energy::EnergyModel& em);
penergy_t FastHairpin(int st, int en);
penergy_t FastTwoLoop(int ost, int oen, int ist, int ien);
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
