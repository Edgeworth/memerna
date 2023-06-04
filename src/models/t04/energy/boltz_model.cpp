// Copyright 2022 Eliot Courtney.
#include "models/t04/energy/boltz_model.h"

#include <memory>
#include <unordered_map>

#include "api/energy/energy_cfg.h"
#include "models/t04/energy/parse.h"

namespace mrna::md::t04 {

BoltzModel::BoltzModel(const Model::Ptr& em) : em_(*em) {
  // Force this to be false to not include bulge states for the partition
  // function.
  em_.cfg.bulge_states = false;
  T04LoadBoltz(*this);
}

}  // namespace mrna::md::t04
