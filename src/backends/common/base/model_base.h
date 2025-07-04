// Copyright 2016 E.
#ifndef BACKENDS_COMMON_BASE_MODEL_BASE_H_
#define BACKENDS_COMMON_BASE_MODEL_BASE_H_

#include <cassert>
#include <cmath>
#include <string>
#include <unordered_map>

#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

using mrna::erg::EnergyResult;

class ModelBase {
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
  std::unordered_map<std::string, Energy> hairpin;
  // Multiloop model:
  Energy multiloop_a = {};
  Energy multiloop_b = {};
  Energy multiloop_c = {};
  // Dangles:
  Energy dangle5[4][4][4] = {};
  Energy dangle3[4][4][4] = {};
  // Coaxial stacking:
  Energy coax_mismatch_non_contiguous = {};
  Energy coax_mismatch_wc_bonus = {};
  Energy coax_mismatch_gu_bonus = {};
  Energy au_penalty = {};
  Energy gu_penalty = {};

  [[nodiscard]] const erg::EnergyCfg& cfg() const { return cfg_; }

  void SetEnergyCfg(const erg::EnergyCfg& cfg) { cfg_ = cfg; }

  [[nodiscard]] constexpr bool CanPair(const Primary& r, int st, int en) const {
    if (cfg_.lonely_pairs == erg::EnergyCfg::LonelyPairs::ON)
      return IsPair(r[st], r[en]) && (en - st - 1 >= HAIRPIN_MIN_SZ);
    return IsPair(r[st], r[en]) && (en - st - 1 >= HAIRPIN_MIN_SZ) &&
        ((en - st - 3 >= HAIRPIN_MIN_SZ && IsPair(r[st + 1], r[en - 1])) ||
            (st > 0 && en < static_cast<int>(r.size() - 1) && IsPair(r[st - 1], r[en + 1])));
  }

  [[nodiscard]] constexpr Energy HairpinInitiation(int n) const {
    assert(n >= 3);
    if (n < INITIATION_CACHE_SZ) return hairpin_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(9) + 1.75 * R * T * ln(n / 9)  -- we use 30 here though to match
    // RNAstructure.
    return hairpin_init[30] + E(1.75 * R * T * log(n / 30.0));
  }

  [[nodiscard]] constexpr Energy BulgeInitiation(int n) const {
    assert(n >= 1);
    if (n < INITIATION_CACHE_SZ) return bulge_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(6) + 1.75 * R * T * ln(n / 6) -- we use 30 here though to match RNAstructure.
    return bulge_init[30] + E(1.75 * R * T * log(n / 30.0));
  }

  [[nodiscard]] constexpr Energy InternalLoopInitiation(int n) const {
    assert(n >= 4);
    if (n < INITIATION_CACHE_SZ) return internal_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    // Formula: G_init(6) + 1.08 * ln(n / 6) -- we use 30 here though to match RNAstructure.
    return internal_init[30] + E(1.08 * log(n / 30.0));
  }

  [[nodiscard]] constexpr Energy MultiloopInitiation(int num_branches, int num_unpaired) const {
    return multiloop_a + num_branches * multiloop_b + num_unpaired * multiloop_c;
  }

  [[nodiscard]] constexpr Energy AuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    if (IsAuPair(stb, enb)) return au_penalty;
    if (IsGuPair(stb, enb)) return gu_penalty;
    return ZERO_E;
  }

  [[nodiscard]] constexpr Energy InternalLoopAuGuPenalty(Base stb, Base enb) const {
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
  [[nodiscard]] Energy MismatchCoaxial(
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

 protected:
  // This is protected to prevent construction on the stack, since this structure is large.
  ModelBase() = default;

  static void Verify(const Primary& r, const Secondary& s, const Ctds* given_ctd) {
    verify(r.size() == s.size(), "sequence and secondary structure must be the same length");
    if (given_ctd)
      verify(given_ctd->size() == r.size(), "given CTDs must be the same length as the seq");
  }

 private:
  mrna::erg::EnergyCfg cfg_ = {};
};

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_MODEL_BASE_H_
