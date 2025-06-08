// Copyright 2025 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_ENERGY_BOLTZ_MODEL_H_
#define BACKENDS_BASEOPT_ENERGY_BOLTZ_MODEL_H_

#include <cassert>

#include "backends/baseopt/energy/model.h"
#include "backends/common/base/boltz_model_base.h"
#include "backends/common/model_mixin.h"

namespace mrna::md::base::opt {

class BoltzModel : public BoltzModelBase<Model>, public ModelMixin<BoltzModel> {
 public:
  BoltzModel() = delete;

  // ModelMixin:
  static BoltzModel::Ptr Create(const Model::Ptr& m) { return BoltzModel::Ptr(new BoltzModel(m)); }

 private:
  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const Model::Ptr& m) : BoltzModelBase(m) {}
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_ENERGY_BOLTZ_MODEL_H_
