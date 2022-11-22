#ifndef COMPUTE_PARTITION_T04_PARTITION_H_
#define COMPUTE_PARTITION_T04_PARTITION_H_

#include <cassert>
#include <tuple>

#include "compute/boltz_dp.h"
#include "compute/energy/energy.h"
#include "compute/energy/t04/boltz_model.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::part::t04 {

std::tuple<BoltzDpArray, BoltzExtArray> PartitionSlowest(
    const Primary& r, const energy::t04::ModelPtr& em);

std::tuple<BoltzDpArray, BoltzExtArray> PartitionFastest(
    const Primary& r, const energy::t04::BoltzModelPtr& bem);

BoltzExtArray Exterior(const Primary& r, const energy::t04::Model& em, const BoltzDpArray& dp);

}  // namespace mrna::part::t04

#endif  // COMPUTE_PARTITION_T04_PARTITION_H_
