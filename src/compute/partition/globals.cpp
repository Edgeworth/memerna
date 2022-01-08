// Copyright 2016 Eliot Courtney.

#include "compute/partition/globals.h"

namespace mrna::partition {

namespace internal {

Array3D<BoltzEnergy, PT_SIZE> gpt;
Array2D<BoltzEnergy, PTEXT_SIZE> gptext;

}  // namespace internal

void SetPartitionGlobalState(const Primary& r) {
  // 0.0 is zero'd memory. TODO: Is this true for mpfr?
  internal::gpt = Array3D<BoltzEnergy, PT_SIZE>(r.size() + 1, 0);
  internal::gptext = Array2D<BoltzEnergy, PTEXT_SIZE>(r.size() + 1, 0);
}

}  // namespace mrna::partition
