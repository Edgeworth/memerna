// Copyright 2016 E.
#include "model/globals.h"

#include "compute/energy/globals.h"
#include "compute/partition/globals.h"

namespace mrna::partition {

namespace internal {

array3d_t<penergy_t, PT_SIZE> gpt;
array2d_t<penergy_t, PTEXT_SIZE> gptext;
precomp_t gppc;

}  // namespace internal

void SetPartitionGlobalState(const primary_t& r, const energy::EnergyModel& em) {
  energy::SetEnergyGlobalState(r, em);
  // 0.0 is zero'd memory.
  internal::gpt = array3d_t<penergy_t, PT_SIZE>(gr.size() + 1, 0);
  internal::gptext = array2d_t<penergy_t, PTEXT_SIZE>(gr.size() + 1, 0);
  internal::gppc = internal::PrecomputeData(r, em);
}

}  // namespace mrna::partition