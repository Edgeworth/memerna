// Copyright 2022 Eliot Courtney.
#include "compute/energy/t04/boltz_model.h"

#include <string>
#include <utility>

#include "compute/energy/common/boltz.h"
#include "compute/energy/energy_cfg.h"

namespace mrna::erg::t04 {

BoltzModel::BoltzModel(const Model::Ptr& em) : T04BoltzMixin(*em) {}

}  // namespace mrna::erg::t04
