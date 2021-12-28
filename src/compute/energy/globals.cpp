// Copyright 2016 Eliot Courtney.
#include "compute/energy/globals.h"

namespace mrna::energy {

EnergyModel gem;
Precomp gpc;

void SetEnergyGlobalState(const Primary& r, const EnergyModel& em) {
  SetGlobalState(r);
  gem = em;
  gpc = PrecomputeData(gr, gem);
}

}  // namespace mrna::energy
