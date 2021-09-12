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
#ifndef MEMERNA_ENERGY_MODEL_H_
#define MEMERNA_ENERGY_MODEL_H_

#include <cmath>
#include <cstdarg>
#include <cstring>

#include "argparse.h"
#include "base.h"
#include "common.h"

namespace memerna {
namespace energy {

class Structure;

class EnergyModel {
public:
  static const int INITIATION_CACHE_SZ = 31;
  // Stacking related:
  // Note that the order of indices is always from 5' to 3'.
  energy_t stack[4][4][4][4];
  // Terminal mismatch:
  energy_t terminal[4][4][4][4];
  // Internal loop related:
  energy_t internal_init[INITIATION_CACHE_SZ];
  energy_t internal_1x1[4][4][4][4][4][4];
  energy_t internal_1x2[4][4][4][4][4][4][4];
  energy_t internal_2x2[4][4][4][4][4][4][4][4];
  energy_t internal_2x3_mismatch[4][4][4][4];
  energy_t internal_other_mismatch[4][4][4][4];
  energy_t internal_asym;
  energy_t internal_augu_penalty;
  // Bulge loop related:
  energy_t bulge_init[INITIATION_CACHE_SZ];
  energy_t bulge_special_c;
  // Hairpin loop related:
  energy_t hairpin_init[INITIATION_CACHE_SZ];
  energy_t hairpin_uu_ga_first_mismatch, hairpin_gg_first_mismatch, hairpin_special_gu_closure,
      hairpin_c3_loop, hairpin_all_c_a, hairpin_all_c_b;
  std::unordered_map<std::string, energy_t> hairpin;
  // Multiloop hack model:
  energy_t multiloop_hack_a, multiloop_hack_b;
  // Dangles:
  energy_t dangle5[4][4][4];
  energy_t dangle3[4][4][4];
  // Coaxial stacking:
  energy_t coax_mismatch_non_contiguous, coax_mismatch_wc_bonus, coax_mismatch_gu_bonus;
  energy_t augu_penalty;

  EnergyModel()
      : stack(), terminal(), internal_init(), internal_1x1(), internal_1x2(), internal_2x2(),
        internal_2x3_mismatch(), internal_other_mismatch(), internal_asym(),
        internal_augu_penalty(), bulge_init(), bulge_special_c(), hairpin_init(),
        hairpin_uu_ga_first_mismatch(), hairpin_gg_first_mismatch(), hairpin_special_gu_closure(),
        hairpin_c3_loop(), hairpin_all_c_a(), hairpin_all_c_b(), hairpin(), multiloop_hack_a(),
        multiloop_hack_b(), dangle5(), dangle3(), coax_mismatch_non_contiguous(),
        coax_mismatch_wc_bonus(), coax_mismatch_gu_bonus(), augu_penalty() {}

  energy_t HairpinInitiation(int n) const {
    assert(n >= 3);
    if (n < INITIATION_CACHE_SZ) return hairpin_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(9) + 1.75 * R * T * ln(n / 9)  -- we use 30 here though to match
    // RNAstructure.
    return energy_t(round(hairpin_init[30] + 10.0 * 1.75 * R * T * log(n / 30.0)));
  }

  energy_t BulgeInitiation(int n) const {
    assert(n >= 1);
    if (n < INITIATION_CACHE_SZ) return bulge_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(6) + 1.75 * R * T * ln(n / 6) -- we use 30 here though to match RNAstructure.
    return energy_t(round(bulge_init[30] + 10.0 * 1.75 * R * T * log(n / 30.0)));
  }

  energy_t InternalLoopInitiation(int n) const {
    assert(n >= 4);
    if (n < INITIATION_CACHE_SZ) return internal_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(6) + 1.08 * ln(n / 6) -- we use 30 here though to match RNAstructure.
    return energy_t(round(internal_init[30] + 10.0 * 1.08 * log(n / 30.0)));
  }

  energy_t MultiloopInitiation(int num_branches) const {
    return multiloop_hack_a + num_branches * multiloop_hack_b;
  }

  energy_t AuGuPenalty(base_t stb, base_t enb) const {
    assert(IsBase(stb) && IsBase(enb));
    return IsAuGu(stb, enb) ? augu_penalty : 0;
  }

  energy_t InternalLoopAuGuPenalty(base_t stb, base_t enb) const {
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
  energy_t MismatchCoaxial(
      base_t five_top, base_t mismatch_top, base_t mismatch_bot, base_t three_bot) const {
    assert(IsBase(five_top) && IsBase(mismatch_top) && IsBase(mismatch_bot) && IsBase(three_bot));
    energy_t coax =
        terminal[five_top][mismatch_top][mismatch_bot][three_bot] + coax_mismatch_non_contiguous;
    if (IsWatsonCrick(mismatch_top, mismatch_bot))
      coax += coax_mismatch_wc_bonus;
    else if (IsGu(mismatch_top, mismatch_bot))
      coax += coax_mismatch_gu_bonus;
    return coax;
  }

  energy_t Hairpin(
      const primary_t& r, int st, int en, std::unique_ptr<Structure>* s = nullptr) const;
  energy_t Bulge(const primary_t& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  energy_t InternalLoop(const primary_t& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  energy_t TwoLoop(const primary_t& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;

  bool IsValid(std::string* reason = nullptr) const;
  uint32_t Checksum() const;
};

typedef std::shared_ptr<EnergyModel> EnergyModelPtr;
}  // namespace energy
}  // namespace memerna

#endif  // MEMERNA_ENERGY_MODEL_H_
