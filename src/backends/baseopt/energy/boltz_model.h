// Copyright 2025 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_ENERGY_BOLTZ_MODEL_H_
#define BACKENDS_BASEOPT_ENERGY_BOLTZ_MODEL_H_

#include <cassert>
#include <memory>

#include "backends/baseopt/energy/model.h"
#include "backends/common/base/boltz_model_base.h"
#include "backends/common/model_mixin.h"

namespace mrna::md::base::opt {

class BoltzModel : public BoltzModelBase<Model>, public ModelMixin<BoltzModel> {
 public:
  BoltzModel() = delete;

  // ModelMixin:
  static BoltzModel::Ptr Create(const Model::Ptr& m) { return BoltzModel::Ptr(new BoltzModel(m)); }

  BoltzEnergy Hairpin(
      const Primary& r, int st, int en, std::unique_ptr<Structure>* s = nullptr) const {
    return m_.Hairpin(r, st, en, s).Boltz();
  }

  BoltzEnergy Bulge(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return m_.Bulge(r, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy InternalLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return m_.InternalLoop(r, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy TwoLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return m_.TwoLoop(r, ost, oen, ist, ien, s).Boltz();
  }

 private:
  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const Model::Ptr& m) : BoltzModelBase(m) {}
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_ENERGY_BOLTZ_MODEL_H_
