// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASE_ENERGY_BOLTZ_MODEL_H_
#define BACKENDS_BASE_ENERGY_BOLTZ_MODEL_H_

#include <cassert>
#include <memory>

#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/boltz_model_base.h"
#include "backends/common/model_mixin.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

class BoltzModel : public BoltzModelBase<Model>, public ModelMixin<BoltzModel> {
 public:
  BoltzModel() = delete;

  // ModelMixin:
  static BoltzModel::Ptr Create(const Model::Ptr& m) { return BoltzModel::Ptr(new BoltzModel(m)); }

  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  BoltzEnergy Hairpin(const Primary& r, const PseudofreeCfgT& pf, int st, int en,
      std::unique_ptr<Structure>* s = nullptr) const {
    return Boltz(m_.Hairpin<CFG>(r, pf, st, en, s));
  }

  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  BoltzEnergy Bulge(const Primary& r, const PseudofreeCfgT& pf, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return Boltz(m_.Bulge<CFG>(r, pf, ost, oen, ist, ien, s));
  }

  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  BoltzEnergy InternalLoop(const Primary& r, const PseudofreeCfgT& pf, int ost, int oen, int ist,
      int ien, std::unique_ptr<Structure>* s = nullptr) const {
    return Boltz(m_.InternalLoop<CFG>(r, pf, ost, oen, ist, ien, s));
  }

  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  BoltzEnergy TwoLoop(const Primary& r, const PseudofreeCfgT& pf, int ost, int oen, int ist,
      int ien, std::unique_ptr<Structure>* s = nullptr) const {
    return Boltz(m_.TwoLoop<CFG>(r, pf, ost, oen, ist, ien, s));
  }

 private:
  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const Model::Ptr& m) : BoltzModelBase(m) {}
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_ENERGY_BOLTZ_MODEL_H_
