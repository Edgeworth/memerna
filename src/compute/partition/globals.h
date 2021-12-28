// Copyright 2016 E.
#ifndef COMPUTE_PARTITION_GLOBALS_H_
#define COMPUTE_PARTITION_GLOBALS_H_

#include "compute/energy/model.h"
#include "compute/partition/partition.h"

namespace mrna::partition {

namespace internal {

extern Array3D<PEnergy, PT_SIZE> gpt;
extern Array2D<PEnergy, PTEXT_SIZE> gptext;
extern Precomp gppc;

}  // namespace internal

void SetPartitionGlobalState(const Primary& r, const energy::EnergyModel& em);

}  // namespace mrna::partition

#endif  // COMPUTE_PARTITION_GLOBALS_H_
