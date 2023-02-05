// Copyright 2023 Eliot Courtney.
#ifndef COMPUTE_ENERGY_T22_BOLTZ_MODEL_H_
#define COMPUTE_ENERGY_T22_BOLTZ_MODEL_H_

#include "model/energy.h"
#include "models/common/model.h"
#include "models/t04/energy/boltz_mixin.h"
#include "models/t22/energy/model.h"

namespace mrna::md::t22 {

class BoltzModel : public ModelMixin<BoltzModel>, public t04::T04BoltzMixin<Model> {
 public:
  BoltzEnergy penultimate_stack[4][4][4][4]{};

  static BoltzModel::Ptr Create(const Model::Ptr& em) {
    return BoltzModel::Ptr(new BoltzModel(em));
  }

 private:
  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const Model::Ptr& em);
};

}  // namespace mrna::md::t22

#endif  // COMPUTE_ENERGY_T22_BOLTZ_MODEL_H_
