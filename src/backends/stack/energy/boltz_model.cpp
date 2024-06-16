// Copyright 2023 E.
#include "backends/stack/energy/boltz_model.h"

#include "backends/common/boltz.h"

namespace mrna::md::stack {

BoltzModel::BoltzModel(const Model::Ptr& m) : BoltzModelBase(m) {
  FILL_BOLTZ((*this), penultimate_stack);
}

}  // namespace mrna::md::stack
