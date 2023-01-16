// Copyright 2023 Eliot Courtney.
#include "compute/energy/t22/boltz_model.h"

#include <string>
#include <utility>

#include "compute/energy/common/boltz.h"
#include "compute/energy/energy_cfg.h"

namespace mrna::erg::t22 {

BoltzModel::BoltzModel(const Model::Ptr& em_unshadow) : T04BoltzMixin(*em_unshadow) {
  FILL_BOLTZ(penultimate_stack);
}

}  // namespace mrna::erg::t22
