// Copyright 2023 Eliot Courtney.
#ifndef COMPUTE_ENERGY_T22_BOLTZ_MODEL_H_
#define COMPUTE_ENERGY_T22_BOLTZ_MODEL_H_

#include "compute/energy/common/model.h"
#include "compute/energy/common/t04like/boltz_mixin.h"
#include "compute/energy/t22/model.h"
#include "model/energy.h"

namespace mrna::erg::t22 {

class BoltzModel : public ModelMixin<BoltzModel>, public T04BoltzMixin<Model> {
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

}  // namespace mrna::erg::t22

#endif  // COMPUTE_ENERGY_T22_BOLTZ_MODEL_H_
