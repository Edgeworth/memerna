// Copyright 2023 E.
#include "models/t22/energy/boltz_model.h"

#include <memory>
#include <unordered_map>

#include "api/energy/energy_cfg.h"
#include "models/common/boltz.h"
#include "models/t04/energy/parse.h"

namespace mrna::md::t22 {

BoltzModel::BoltzModel(const Model::Ptr& em) : em_(*em) {
  // Force this to be false to not include bulge states for the partition
  // function.
  em_.cfg.bulge_states = false;
  t04::T04LoadBoltz(*this);

  FILL_BOLTZ((*this), penultimate_stack);
}

}  // namespace mrna::md::t22
