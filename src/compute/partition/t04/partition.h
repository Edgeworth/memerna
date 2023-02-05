// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_PARTITION_T04_PARTITION_H_
#define COMPUTE_PARTITION_T04_PARTITION_H_

#include <cassert>
#include <tuple>

#include "compute/energy/energy.h"
#include "compute/energy/t04/boltz_model.h"
#include "compute/partition/t04/dp.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::part::t04 {

std::tuple<BoltzDpArray, BoltzExtArray> PartitionSlowest(
    const Primary& r, const erg::t04::Model::Ptr& em);

std::tuple<BoltzDpArray, BoltzExtArray> PartitionFastest(
    const Primary& r, const erg::t04::BoltzModel::Ptr& bem);

BoltzExtArray PartitionExterior(
    const Primary& r, const erg::t04::Model& em, const BoltzDpArray& dp);

}  // namespace mrna::part::t04

#endif  // COMPUTE_PARTITION_T04_PARTITION_H_
