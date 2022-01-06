// Copyright 2016 Eliot Courtney.
#include "compute/mfe/globals.h"

#include "compute/constants.h"
#include "compute/dp.h"
#include "compute/energy/globals.h"

namespace mrna::mfe {

namespace internal {

std::vector<int> gp;
std::vector<Ctd> gctd;
std::string grep;
Energy genergy;
Array3D<Energy, DP_SIZE> gdp;
Array2D<Energy, EXT_SIZE> gext;

}  // namespace internal

void SetMfeGlobalState(const Primary& r, const energy::EnergyModel& em) {
  energy::SetEnergyGlobalState(r, em);
  internal::gp.resize(r.size());
  internal::gctd.resize(r.size());
  internal::genergy = MAX_E;
  std::fill(internal::gp.begin(), internal::gp.end(), -1);
  std::fill(internal::gctd.begin(), internal::gctd.end(), CTD_NA);
  internal::gdp = Array3D<Energy, DP_SIZE>(r.size() + 1);
  // TODO: Move this
  internal::gext = Array2D<Energy, EXT_SIZE>(r.size() + 1);
}

}  // namespace mrna::mfe
