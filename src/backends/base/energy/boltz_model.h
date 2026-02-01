// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASE_ENERGY_BOLTZ_MODEL_H_
#define BACKENDS_BASE_ENERGY_BOLTZ_MODEL_H_

#include <cassert>
#include <memory>

#include "backends/base/energy/model.h"
#include "backends/common/base/boltz_model_base.h"
#include "backends/common/model_mixin.h"

namespace mrna::md::base {

class BoltzModel : public BoltzModelBase<Model>, public ModelMixin<BoltzModel> {
 public:
  BoltzModel() = delete;

  // ModelMixin:
  static BoltzModel::Ptr Create(const Model::Ptr& m) { return BoltzModel::Ptr(new BoltzModel(m)); }

  BoltzEnergy Hairpin(const Primary& r, const erg::PseudofreeCfg& pf, int st, int en,
      std::unique_ptr<Structure>* s = nullptr) const {
    return m_.Hairpin(r, pf, st, en, s).Boltz();
  }

  BoltzEnergy Bulge(const Primary& r, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, int ost,
      int oen, int ist, int ien, std::unique_ptr<Structure>* s = nullptr) const {
    return m_.Bulge(r, cfg, pf, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy InternalLoop(const Primary& r, const erg::PseudofreeCfg& pf, int ost, int oen,
      int ist, int ien, std::unique_ptr<Structure>* s = nullptr) const {
    return m_.InternalLoop(r, pf, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy TwoLoop(const Primary& r, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, int ost,
      int oen, int ist, int ien, std::unique_ptr<Structure>* s = nullptr) const {
    return m_.TwoLoop(r, cfg, pf, ost, oen, ist, ien, s).Boltz();
  }

 private:
  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const Model::Ptr& m) : BoltzModelBase(m) {}
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_ENERGY_BOLTZ_MODEL_H_
