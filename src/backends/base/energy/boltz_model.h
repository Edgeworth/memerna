// Copyright 2022 E.
#ifndef BACKENDS_BASE_ENERGY_BOLTZ_MODEL_H_
#define BACKENDS_BASE_ENERGY_BOLTZ_MODEL_H_

#include <cassert>

#include "backends/base/energy/model.h"
#include "backends/common/base/boltz_model_base.h"
#include "backends/common/base/boltz_pseudofree_model.h"
#include "backends/common/model_mixin.h"

namespace mrna::md::base {

class BoltzModel : public BoltzModelBase<Model>, public ModelMixin<BoltzModel> {
 public:
  BoltzPseudofreeModel pf;

  BoltzModel() = delete;

  // ModelMixin:
  static BoltzModel::Ptr Create(const Model::Ptr& m) { return BoltzModel::Ptr(new BoltzModel(m)); }

 private:
  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModel(const Model::Ptr& m) : BoltzModelBase(m) {
    pf.Load(m->pf.paired, m->pf.unpaired);
  }
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_ENERGY_BOLTZ_MODEL_H_
