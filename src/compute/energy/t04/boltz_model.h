// Copyright 2022 E.
#ifndef COMPUTE_ENERGY_BOLTZMANN_MODEL_H_
#define COMPUTE_ENERGY_BOLTZMANN_MODEL_H_

#include <cassert>
#include <memory>
#include <string>
#include <unordered_map>

#include "compute/energy/t04/model.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::erg::t04 {

class BoltzModel;
using BoltzModelPtr = std::shared_ptr<BoltzModel>;

class BoltzModel {
 public:
  BoltzEnergy stack[4][4][4][4]{};
  BoltzEnergy terminal[4][4][4][4]{};
  BoltzEnergy internal_init[Model::INITIATION_CACHE_SZ]{};
  BoltzEnergy internal_1x1[4][4][4][4][4][4]{};
  BoltzEnergy internal_1x2[4][4][4][4][4][4][4]{};
  BoltzEnergy internal_2x2[4][4][4][4][4][4][4][4]{};
  BoltzEnergy internal_2x3_mismatch[4][4][4][4]{};
  BoltzEnergy internal_other_mismatch[4][4][4][4]{};
  BoltzEnergy internal_asym{};
  BoltzEnergy internal_augu_penalty{};
  BoltzEnergy bulge_init[Model::INITIATION_CACHE_SZ]{};
  BoltzEnergy bulge_special_c{};
  BoltzEnergy hairpin_init[Model::INITIATION_CACHE_SZ]{};
  BoltzEnergy hairpin_uu_ga_first_mismatch{}, hairpin_gg_first_mismatch{},
      hairpin_special_gu_closure{}, hairpin_c3_loop{}, hairpin_all_c_a{}, hairpin_all_c_b{};
  std::unordered_map<std::string, BoltzEnergy> hairpin;
  BoltzEnergy multiloop_hack_a{}, multiloop_hack_b{};
  BoltzEnergy dangle5[4][4][4]{};
  BoltzEnergy dangle3[4][4][4]{};
  BoltzEnergy coax_mismatch_non_contiguous{}, coax_mismatch_wc_bonus{}, coax_mismatch_gu_bonus{};
  BoltzEnergy augu_penalty{};

  BoltzModel() = delete;

  static BoltzModelPtr Create(const ModelPtr& em) { return BoltzModelPtr(new BoltzModel(em)); }

  const ModelPtr& em() const { return em_; }

  BoltzEnergy InternalLoopAuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    return IsAuGuPair(stb, enb) ? internal_augu_penalty : 1.0;
  }

  BoltzEnergy AuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    return IsAuGuPair(stb, enb) ? augu_penalty : 1.0;
  }

  BoltzEnergy MismatchCoaxial(
      Base five_top, Base mismatch_top, Base mismatch_bot, Base three_bot) const {
    assert(IsBase(five_top) && IsBase(mismatch_top) && IsBase(mismatch_bot) && IsBase(three_bot));
    BoltzEnergy coax =
        terminal[five_top][mismatch_top][mismatch_bot][three_bot] * coax_mismatch_non_contiguous;
    if (IsWcPair(mismatch_top, mismatch_bot))
      coax *= coax_mismatch_wc_bonus;
    else if (IsGuPair(mismatch_top, mismatch_bot))
      coax *= coax_mismatch_gu_bonus;
    return coax;
  }

  // TODO(2): Implement versions of Bulge, InternalLoop, TwoLoop, Hairpin with boltzmann baked in.
  BoltzEnergy Hairpin(
      const Primary& r, int st, int en, std::unique_ptr<Structure>* s = nullptr) const {
    return em()->Hairpin(r, st, en, s).Boltz();
  }

  BoltzEnergy Bulge(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return em()->Bulge(r, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy InternalLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return em()->InternalLoop(r, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy TwoLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return em()->TwoLoop(r, ost, oen, ist, ien, s).Boltz();
  }

 private:
  ModelPtr em_;

  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const ModelPtr& em);
};

}  // namespace mrna::erg::t04

#endif  // COMPUTE_ENERGY_BOLTZMANN_MODEL_H_
