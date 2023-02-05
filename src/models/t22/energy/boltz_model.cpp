// Copyright 2023 Eliot Courtney.
#include "models/t22/energy/boltz_model.h"

#include <memory>
#include <unordered_map>

#include "models/common/boltz.h"

namespace mrna::md::t22::erg {

BoltzModel::BoltzModel(const Model::Ptr& em_unshadow) : T04BoltzMixin(*em_unshadow) {
  FILL_BOLTZ(penultimate_stack);
}

}  // namespace mrna::md::t22::erg
