// Copyright 2016 Eliot Courtney.
#ifndef PARTITION_PARTITION_GLOBALS_H_
#define PARTITION_PARTITION_GLOBALS_H_

#include "common.h"
#include "energy/energy_model.h"
#include "partition/partition.h"

namespace mrna::partition {

namespace internal {

extern array3d_t<penergy_t, PT_SIZE> gpt;
extern array2d_t<penergy_t, PTEXT_SIZE> gptext;
extern precomp_t gppc;

}  // namespace internal

void SetPartitionGlobalState(const primary_t& r, const energy::EnergyModel& em);

}  // namespace mrna::partition

#endif  // PARTITION_PARTITION_GLOBALS_H_
