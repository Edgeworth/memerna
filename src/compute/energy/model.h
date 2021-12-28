// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_MODEL_H_
#define COMPUTE_ENERGY_MODEL_H_

#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_map>

#include "model/base.h"
#include "model/structure.h"
#include "util/argparse.h"

namespace mrna::energy {

class Structure;

class EnergyModel {
 public:
  static const int INITIATION_CACHE_SZ = 31;
  // Stacking related:
  // Note that the order of indices is always from 5' to 3'.
  Energy stack[4][4][4][4];
  // Terminal mismatch:
  Energy terminal[4][4][4][4];
  // Internal loop related:
  Energy internal_init[INITIATION_CACHE_SZ];
  Energy internal_1x1[4][4][4][4][4][4];
  Energy internal_1x2[4][4][4][4][4][4][4];
  Energy internal_2x2[4][4][4][4][4][4][4][4];
  Energy internal_2x3_mismatch[4][4][4][4];
  Energy internal_other_mismatch[4][4][4][4];
  Energy internal_asym;
  Energy internal_augu_penalty;
  // Bulge loop related:
  Energy bulge_init[INITIATION_CACHE_SZ];
  Energy bulge_special_c;
  // Hairpin loop related:
  Energy hairpin_init[INITIATION_CACHE_SZ];
  Energy hairpin_uu_ga_first_mismatch, hairpin_gg_first_mismatch, hairpin_special_gu_closure,
      hairpin_c3_loop, hairpin_all_c_a, hairpin_all_c_b;
  std::unordered_map<std::string, Energy> hairpin;
  // Multiloop hack model:
  Energy multiloop_hack_a, multiloop_hack_b;
  // Dangles:
  Energy dangle5[4][4][4];
  Energy dangle3[4][4][4];
  // Coaxial stacking:
  Energy coax_mismatch_non_contiguous, coax_mismatch_wc_bonus, coax_mismatch_gu_bonus;
  Energy augu_penalty;

  EnergyModel()
      : stack(), terminal(), internal_init(), internal_1x1(), internal_1x2(), internal_2x2(),
        internal_2x3_mismatch(), internal_other_mismatch(), internal_asym(),
        internal_augu_penalty(), bulge_init(), bulge_special_c(), hairpin_init(),
        hairpin_uu_ga_first_mismatch(), hairpin_gg_first_mismatch(), hairpin_special_gu_closure(),
        hairpin_c3_loop(), hairpin_all_c_a(), hairpin_all_c_b(), hairpin(), multiloop_hack_a(),
        multiloop_hack_b(), dangle5(), dangle3(), coax_mismatch_non_contiguous(),
        coax_mismatch_wc_bonus(), coax_mismatch_gu_bonus(), augu_penalty() {}

  Energy HairpinInitiation(int n) const {
    assert(n >= 3);
    if (n < INITIATION_CACHE_SZ) return hairpin_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(9) + 1.75 * R * T * ln(n / 9)  -- we use 30 here though to match
    // RNAstructure.
    return Energy(round(hairpin_init[30] + 10.0 * 1.75 * R * T * log(n / 30.0)));
  }

  Energy BulgeInitiation(int n) const {
    assert(n >= 1);
    if (n < INITIATION_CACHE_SZ) return bulge_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(6) + 1.75 * R * T * ln(n / 6) -- we use 30 here though to match RNAstructure.
    return Energy(round(bulge_init[30] + 10.0 * 1.75 * R * T * log(n / 30.0)));
  }

  Energy InternalLoopInitiation(int n) const {
    assert(n >= 4);
    if (n < INITIATION_CACHE_SZ) return internal_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(6) + 1.08 * ln(n / 6) -- we use 30 here though to match RNAstructure.
    return Energy(round(internal_init[30] + 10.0 * 1.08 * log(n / 30.0)));
  }

  Energy MultiloopInitiation(int num_branches) const {
    return multiloop_hack_a + num_branches * multiloop_hack_b;
  }

  Energy AuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    return IsAuGu(stb, enb) ? augu_penalty : 0;
  }

  Energy InternalLoopAuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    return IsAuGu(stb, enb) ? internal_augu_penalty : 0;
  }

  // We use the normal terminal mismatch parameters for the mismatch that is on the continuous part
  // of the RNA. The stacking for the non-continuous part is set to be an arbitrary given number.
  // There are two possible orientations, since the base involved in the terminal mismatch could
  // come from either side.
  // ... _ _ _ _ ...
  // ...|_|  _|_|...
  //      | |
  // Rules for mismatch mediated coaxial stacking:
  //    1. A terminal mismatch is formed around the branch being straddled.
  //    2. An arbitrary bonus is added.
  //    2. An arbitrary bonus is added if the mismatch is Watson-Crick or GU.
  Energy MismatchCoaxial(
      Base five_top, Base mismatch_top, Base mismatch_bot, Base three_bot) const {
    assert(IsBase(five_top) && IsBase(mismatch_top) && IsBase(mismatch_bot) && IsBase(three_bot));
    Energy coax =
        terminal[five_top][mismatch_top][mismatch_bot][three_bot] + coax_mismatch_non_contiguous;
    if (IsWatsonCrick(mismatch_top, mismatch_bot))
      coax += coax_mismatch_wc_bonus;
    else if (IsGu(mismatch_top, mismatch_bot))
      coax += coax_mismatch_gu_bonus;
    return coax;
  }

  Energy Hairpin(const Primary& r, int st, int en, std::unique_ptr<Structure>* s = nullptr) const;
  Energy Bulge(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy InternalLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy TwoLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;

  bool IsValid(std::string* reason = nullptr) const;
  uint32_t Checksum() const;
};

typedef std::shared_ptr<EnergyModel> EnergyModelPtr;

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_MODEL_H_
