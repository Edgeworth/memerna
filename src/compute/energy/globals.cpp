// Copyright 2016 E.
#include "compute/energy/globals.h"

namespace mrna::energy {

EnergyModel gem;
Precomp gpc;

void SetEnergyGlobalState(const Primary& r, const EnergyModel& em) {
  gem = em;
  gpc = PrecomputeData(r, gem);
}

}  // namespace mrna::energy
