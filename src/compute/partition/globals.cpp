// Copyright 2016 Eliot Courtney.

#include "compute/energy/globals.h"

#include "compute/partition/globals.h"

namespace mrna::partition {

namespace internal {

Array3D<PEnergy, PT_SIZE> gpt;
Array2D<PEnergy, PTEXT_SIZE> gptext;
Precomp gppc;

}  // namespace internal

void SetPartitionGlobalState(const Primary& r, const energy::EnergyModel& em) {
  energy::SetEnergyGlobalState(r, em);
  // 0.0 is zero'd memory.
  internal::gpt = Array3D<PEnergy, PT_SIZE>(gr.size() + 1, 0);
  internal::gptext = Array2D<PEnergy, PTEXT_SIZE>(gr.size() + 1, 0);
  internal::gppc = internal::PrecomputeData(r, em);
}

}  // namespace mrna::partition
