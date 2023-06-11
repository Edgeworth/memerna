// Copyright 2022 Eliot Courtney.
#ifndef MODELS_T22_ENERGY_MODEL_H_
#define MODELS_T22_ENERGY_MODEL_H_

#include <fmt/core.h>

#include <cassert>
#include <cmath>
#include <deque>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>

#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "model/structure.h"
#include "models/common/model.h"
#include "models/t04/energy/model.h"
#include "models/t04/energy/parse.h"

namespace mrna::md::t22 {

using mrna::erg::EnergyResult;

class Model : public ModelMixin<Model> {
 public:
  static constexpr int INITIATION_CACHE_SZ = t04::Model::INITIATION_CACHE_SZ;
  static constexpr double RAND_MIN_ENERGY = t04::Model::RAND_MIN_ENERGY;
  static constexpr double RAND_MAX_ENERGY = t04::Model::RAND_MAX_ENERGY;
  static constexpr int RAND_MAX_HAIRPIN_SZ = t04::Model::RAND_MAX_HAIRPIN_SZ;
  static constexpr int RAND_MAX_NUM_HAIRPIN = t04::Model::RAND_MAX_NUM_HAIRPIN;

  // From T04:
  Energy stack[4][4][4][4] = {};
  Energy terminal[4][4][4][4] = {};
  Energy internal_init[INITIATION_CACHE_SZ] = {};
  Energy internal_1x1[4][4][4][4][4][4] = {};
  Energy internal_1x2[4][4][4][4][4][4][4] = {};
  Energy internal_2x2[4][4][4][4][4][4][4][4] = {};
  Energy internal_2x3_mismatch[4][4][4][4] = {};
  Energy internal_other_mismatch[4][4][4][4] = {};
  Energy internal_asym = {};
  Energy internal_au_penalty = {};
  Energy internal_gu_penalty = {};
  Energy bulge_init[INITIATION_CACHE_SZ] = {};
  Energy bulge_special_c = {};
  Energy hairpin_init[INITIATION_CACHE_SZ] = {};
  Energy hairpin_uu_ga_first_mismatch = {};
  Energy hairpin_gg_first_mismatch = {};
  Energy hairpin_special_gu_closure = {};
  Energy hairpin_c3_loop = {};
  Energy hairpin_all_c_a = {};
  Energy hairpin_all_c_b = {};
  std::unordered_map<std::string, Energy> hairpin = {};
  Energy multiloop_hack_a = {};
  Energy multiloop_hack_b = {};
  Energy dangle5[4][4][4] = {};
  Energy dangle3[4][4][4] = {};
  Energy coax_mismatch_non_contiguous = {};
  Energy coax_mismatch_wc_bonus = {};
  Energy coax_mismatch_gu_bonus = {};
  Energy au_penalty = {};
  Energy gu_penalty = {};
  // Specific to T22:
  Energy penultimate_stack[4][4][4][4] = {};

  // Pseudofree energies. Ignored if empty.
  std::vector<Energy> pf_paired;
  std::vector<Energy> pf_unpaired;
  // Exclusive cumulative sum of size N+1
  std::vector<Energy> pf_unpaired_cum;

  mrna::erg::EnergyCfg cfg = {};

  [[nodiscard]] inline constexpr bool CanPair(const Primary& r, int st, int en) const {
    if (cfg.lonely_pairs == erg::EnergyCfg::LonelyPairs::ON)
      return IsPair(r[st], r[en]) && (en - st - 1 >= HAIRPIN_MIN_SZ);
    return IsPair(r[st], r[en]) && (en - st - 1 >= HAIRPIN_MIN_SZ) &&
        ((en - st - 3 >= HAIRPIN_MIN_SZ && IsPair(r[st + 1], r[en - 1])) ||
            (st > 0 && en < static_cast<int>(r.size() - 1) && IsPair(r[st - 1], r[en + 1])));
  }

  [[nodiscard]] constexpr Energy HairpinInitiation(int n) const {
    assert(n >= 3);
    if (n < INITIATION_CACHE_SZ) return hairpin_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    return hairpin_init[30] + E(1.75 * R * T * log(n / 30.0));
  }

  [[nodiscard]] constexpr Energy BulgeInitiation(int n) const {
    assert(n >= 1);
    if (n < INITIATION_CACHE_SZ) return bulge_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    return bulge_init[30] + E(1.75 * R * T * log(n / 30.0));
  }

  [[nodiscard]] constexpr Energy InternalLoopInitiation(int n) const {
    assert(n >= 4);
    if (n < INITIATION_CACHE_SZ) return internal_init[n];
    static_assert(INITIATION_CACHE_SZ > 30, "Need initiation values for up to 30.");
    return internal_init[30] + E(1.08 * log(n / 30.0));
  }

  [[nodiscard]] constexpr Energy MultiloopInitiation(int num_branches) const {
    return multiloop_hack_a + num_branches * multiloop_hack_b;
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

  Energy Hairpin(const Primary& r, int st, int en, std::unique_ptr<Structure>* s = nullptr) const;
  Energy Bulge(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy InternalLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy TwoLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy MultiloopEnergy(const Primary& r, const Secondary& s, int st, int en,
      std::deque<int>* branches, bool use_given_ctds, Ctds* ctd,
      std::unique_ptr<Structure>* sstruc = nullptr) const;

  [[nodiscard]] constexpr Energy PfUnpaired(int n) const {
    if (pf_unpaired.empty()) return ZERO_E;
    return pf_unpaired[n];
  }

  // Inclusive range, unlike pf_unpaired_cum directly.
  [[nodiscard]] constexpr Energy PfUnpairedCum(int st, int en) const {
    if (pf_unpaired.empty()) return ZERO_E;
    return pf_unpaired_cum[en + 1] - pf_unpaired_cum[st];
  }

  [[nodiscard]] constexpr Energy PfPaired(int st, int en) const {
    if (pf_paired.empty()) return ZERO_E;
    return pf_paired[st] + pf_paired[en];
  }

  // Computes the penalty for a stack of the given length, ending at (ist, ien).
  // Handles bulge loops.
  [[nodiscard]] constexpr Energy StackPenalty(const Primary& r, const Secondary& s, int ost,
      int oen, int ist, int ien, std::unique_ptr<Structure>* struc = nullptr) const;

  // ModelMixin:
  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
      int en, bool build_structure = false) const;
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      bool build_structure = false) const;

  bool IsValid(std::string* reason = nullptr) const { return t04::T04IsValid(*this, reason); }

  void LoadPseudofreeEnergy(std::vector<Energy> paired, std::vector<Energy> unpaired);

 protected:
  void LoadFromModelPath(const std::string& path);
  void LoadRandom(std::mt19937& eng);

 private:
  friend class ModelMixin<Model>;

  void VerifyLengths(const Primary& r, const Secondary& s, const Ctds* given_ctd) const;

  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;

  Energy SubEnergyInternal(const Primary& r, const Secondary& s, int st, int en, int stack_st,
      int stack_en, bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const;
};

}  // namespace mrna::md::t22

#endif  // MODELS_T22_ENERGY_MODEL_H_
