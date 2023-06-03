// Copyright 2023 Eliot Courtney.
#ifndef MODELS_T22_ENERGY_BOLTZ_MODEL_H_
#define MODELS_T22_ENERGY_BOLTZ_MODEL_H_

#include "model/energy.h"
#include "models/common/model.h"
#include "models/t22/energy/model.h"

namespace mrna::md::t22 {

class BoltzModel : public ModelMixin<BoltzModel> {
 public:
  // From T04:
  BoltzEnergy stack[4][4][4][4]{};
  BoltzEnergy terminal[4][4][4][4]{};
  BoltzEnergy internal_init[Model::INITIATION_CACHE_SZ]{};
  BoltzEnergy internal_1x1[4][4][4][4][4][4]{};
  BoltzEnergy internal_1x2[4][4][4][4][4][4][4]{};
  BoltzEnergy internal_2x2[4][4][4][4][4][4][4][4]{};
  BoltzEnergy internal_2x3_mismatch[4][4][4][4]{};
  BoltzEnergy internal_other_mismatch[4][4][4][4]{};
  BoltzEnergy internal_asym{};
  BoltzEnergy internal_au_penalty{};
  BoltzEnergy internal_gu_penalty{};
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
  BoltzEnergy au_penalty{};
  BoltzEnergy gu_penalty{};
  // Specific to T22:
  BoltzEnergy penultimate_stack[4][4][4][4]{};

  // ModelMixin:
  static BoltzModel::Ptr Create(const Model::Ptr& em) {
    return BoltzModel::Ptr(new BoltzModel(em));
  }

  const Model& em() const { return em_; }

 private:
  Model em_;

  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const Model::Ptr& em);
};

}  // namespace mrna::md::t22

#endif  // MODELS_T22_ENERGY_BOLTZ_MODEL_H_
