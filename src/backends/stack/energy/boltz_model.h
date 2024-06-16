// Copyright 2023 Eliot Courtney.
#ifndef BACKENDS_STACK_ENERGY_BOLTZ_MODEL_H_
#define BACKENDS_STACK_ENERGY_BOLTZ_MODEL_H_

#include "backends/common/base/boltz_model_base.h"
#include "backends/common/model_mixin.h"
#include "backends/stack/energy/model.h"
#include "model/energy.h"

namespace mrna::md::stack {

class BoltzModel : public base::BoltzModelBase<Model>, public ModelMixin<BoltzModel> {
 public:
  BoltzEnergy penultimate_stack[4][4][4][4]{};

  // ModelMixin:
  static BoltzModel::Ptr Create(const Model::Ptr& m) { return BoltzModel::Ptr(new BoltzModel(m)); }

 private:
  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const Model::Ptr& m);
};

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_ENERGY_BOLTZ_MODEL_H_
