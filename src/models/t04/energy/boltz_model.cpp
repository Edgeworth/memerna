// Copyright 2022 E.
#include "models/t04/energy/boltz_model.h"

#include <memory>
#include <unordered_map>

namespace mrna::md::t04::erg {

BoltzModel::BoltzModel(const Model::Ptr& em) : T04BoltzMixin(*em) {}

}  // namespace mrna::md::t04::erg