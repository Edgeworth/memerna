// Copyright 2023 E.
#ifndef COMPUTE_ENERGY_COMMON_T04LIKE_MODEL_MIXIN_H_
#define COMPUTE_ENERGY_COMMON_T04LIKE_MODEL_MIXIN_H_

#include <cassert>
#include <cmath>
#include <random>
#include <string>
#include <unordered_map>

#include "compute/energy/energy_cfg.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::erg {

class T04ModelMixin {
 public:
  static constexpr int INITIATION_CACHE_SZ = 31;
  static constexpr double RAND_MIN_ENERGY = -10.0;
  static constexpr double RAND_MAX_ENERGY = 10.0;
  static constexpr int RAND_MAX_HAIRPIN_SZ = 8;
  static constexpr int RAND_MAX_NUM_HAIRPIN = 50;

  // Stacking related:
  // Note that the order of indices is always from 5' to 3'.
  Energy stack[4][4][4][4] = {};
  // Terminal mismatch:
  Energy terminal[4][4][4][4] = {};
  // Internal loop related:
  Energy internal_init[INITIATION_CACHE_SZ] = {};
  Energy internal_1x1[4][4][4][4][4][4] = {};
  Energy internal_1x2[4][4][4][4][4][4][4] = {};
  Energy internal_2x2[4][4][4][4][4][4][4][4] = {};
  Energy internal_2x3_mismatch[4][4][4][4] = {};
  Energy internal_other_mismatch[4][4][4][4] = {};
  Energy internal_asym = {};
  Energy internal_au_penalty = {};
  Energy internal_gu_penalty = {};
  // Bulge loop related:
  Energy bulge_init[INITIATION_CACHE_SZ] = {};
  Energy bulge_special_c = {};
  // Hairpin loop related:
  Energy hairpin_init[INITIATION_CACHE_SZ] = {};
  Energy hairpin_uu_ga_first_mismatch = {};
  Energy hairpin_gg_first_mismatch = {};
  Energy hairpin_special_gu_closure = {};
  Energy hairpin_c3_loop = {};
  Energy hairpin_all_c_a = {};
  Energy hairpin_all_c_b = {};
  std::unordered_map<std::string, Energy> hairpin = {};
  // Multiloop hack model:
  Energy multiloop_hack_a = {};
  Energy multiloop_hack_b = {};
  // Dangles:
  Energy dangle5[4][4][4] = {};
  Energy dangle3[4][4][4] = {};
  // Coaxial stacking:
  Energy coax_mismatch_non_contiguous = {};
  Energy coax_mismatch_wc_bonus = {};
  Energy coax_mismatch_gu_bonus = {};
  Energy au_penalty = {};
  Energy gu_penalty = {};

  EnergyCfg cfg = {};

  inline bool CanPair(const Primary& r, int st, int en) const {
    if (cfg.lonely_pairs) return IsPair(r[st], r[en]) && (en - st - 1 >= HAIRPIN_MIN_SZ);
    return IsPair(r[st], r[en]) && (en - st - 1 >= HAIRPIN_MIN_SZ) &&
        ((en - st - 3 >= HAIRPIN_MIN_SZ && IsPair(r[st + 1], r[en - 1])) ||
            (st > 0 && en < static_cast<int>(r.size() - 1) && IsPair(r[st - 1], r[en + 1])));
  }

  Energy HairpinInitiation(int n) const {
    assert(n >= 3);
    if (n < INITIATION_CACHE_SZ) return hairpin_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(9) + 1.75 * R * T * ln(n / 9)  -- we use 30 here though to match
    // RNAstructure.
    return hairpin_init[30] + E(1.75 * R * T * log(n / 30.0));
  }

  Energy BulgeInitiation(int n) const {
    assert(n >= 1);
    if (n < INITIATION_CACHE_SZ) return bulge_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(6) + 1.75 * R * T * ln(n / 6) -- we use 30 here though to match RNAstructure.
    return bulge_init[30] + E(1.75 * R * T * log(n / 30.0));
  }

  Energy InternalLoopInitiation(int n) const {
    assert(n >= 4);
    if (n < INITIATION_CACHE_SZ) return internal_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(6) + 1.08 * ln(n / 6) -- we use 30 here though to match RNAstructure.
    return internal_init[30] + E(1.08 * log(n / 30.0));
  }

  Energy MultiloopInitiation(int num_branches) const {
    return multiloop_hack_a + num_branches * multiloop_hack_b;
  }

  Energy AuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    if (IsAuPair(stb, enb)) return au_penalty;
    if (IsGuPair(stb, enb)) return gu_penalty;
    return ZERO_E;
  }

  Energy InternalLoopAuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    if (IsAuPair(stb, enb)) return internal_au_penalty;
    if (IsGuPair(stb, enb)) return internal_gu_penalty;
    return ZERO_E;
  }

  // We use the normal terminal mismatch parameters for the mismatch that is on
  // the continuous part of the RNA. The stacking for the non-continuous part is
  // set to be an arbitrary given number. There are two possible orientations,
  // since the base involved in the terminal mismatch could come from either
  // side.
  // ... _ _ _ _ ...
  // ...|_|  _|_|...
  //      | |
  // Rules for mismatch mediated coaxial stacking:
  // 1. A terminal mismatch is formed around the branch being straddled.
  // 2. An arbitrary bonus is added.
  // 2. An arbitrary bonus is added if the mismatch is Watson-Crick or GU.
  Energy MismatchCoaxial(
      Base five_top, Base mismatch_top, Base mismatch_bot, Base three_bot) const {
    assert(IsBase(five_top) && IsBase(mismatch_top) && IsBase(mismatch_bot) && IsBase(three_bot));
    Energy coax =
        terminal[five_top][mismatch_top][mismatch_bot][three_bot] + coax_mismatch_non_contiguous;
    if (IsWcPair(mismatch_top, mismatch_bot))
      coax += coax_mismatch_wc_bonus;
    else if (IsGuPair(mismatch_top, mismatch_bot))
      coax += coax_mismatch_gu_bonus;
    return coax;
  }

  bool IsValid(std::string* reason = nullptr) const;

 protected:
  void LoadFromDir(const std::string& data_dir);
  void LoadRandom(std::mt19937& eng);

 private:
  void ParseMiscDataFromFile(const std::string& filename);
};

}  // namespace mrna::erg

#endif  // COMPUTE_ENERGY_COMMON_T04LIKE_MODEL_MIXIN_H_
