// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_T04_BOLTZ_MODEL_H_
#define COMPUTE_ENERGY_T04_BOLTZ_MODEL_H_

#include "compute/energy/common/model.h"
#include "compute/energy/t04/boltz_mixin.h"
#include "compute/energy/t04/model.h"

namespace mrna::erg::t04 {

class BoltzModel : public ModelMixin<BoltzModel>, public T04BoltzMixin<Model> {
 public:
  BoltzModel() = delete;

  static BoltzModel::Ptr Create(const Model::Ptr& em) {
    return BoltzModel::Ptr(new BoltzModel(em));
  }

 private:
  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const Model::Ptr& em);
};

}  // namespace mrna::erg::t04

#endif  // COMPUTE_ENERGY_T04_BOLTZ_MODEL_H_
