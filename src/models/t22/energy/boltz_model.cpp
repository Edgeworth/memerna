// Copyright 2023 Eliot Courtney.
#include "compute/energy/t22/boltz_model.h"

#include <memory>
#include <unordered_map>

#include "models/common/boltz.h"

namespace mrna::md::t22::erg {

BoltzModel::BoltzModel(const Model::Ptr& em_unshadow) : T04BoltzMixin(*em_unshadow) {
  FILL_BOLTZ(penultimate_stack);
}

}  // namespace mrna::md::t22::erg
