// Copyright 2016 E.
#include "energy/energy_globals.h"

namespace mrna::energy {

EnergyModel gem;
precomp_t gpc;

void SetEnergyGlobalState(const primary_t& r, const EnergyModel& em) {
  SetGlobalState(r);
  gem = em;
  gpc = PrecomputeData(gr, gem);
}

}  // namespace mrna::energy
