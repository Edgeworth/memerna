// Copyright 2016 E.
#ifndef COMPUTE_PARTITION_GLOBALS_H_
#define COMPUTE_PARTITION_GLOBALS_H_

#include "compute/energy/model.h"
#include "compute/partition/partition.h"

namespace mrna::partition {

namespace internal {

extern Array3D<BoltzEnergy, PT_SIZE> gpt;
extern Array2D<BoltzEnergy, PTEXT_SIZE> gptext;

}  // namespace internal

void SetPartitionGlobalState(const Primary& r);

}  // namespace mrna::partition

#endif  // COMPUTE_PARTITION_GLOBALS_H_
