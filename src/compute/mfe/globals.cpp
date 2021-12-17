// Copyright 2016 E.
#include "compute/mfe/globals.h"

#include "compute/constants.h"
#include "compute/dp.h"
#include "compute/energy/globals.h"
#include "model/globals.h"

namespace mrna::mfe {

namespace internal {

std::vector<int> gp;
std::vector<Ctd> gctd;
std::string grep;
energy_t genergy;
array3d_t<energy_t, DP_SIZE> gdp;
array2d_t<energy_t, EXT_SIZE> gext;

}  // namespace internal

void SetMfeGlobalState(const primary_t& r, const energy::EnergyModel& em) {
  energy::SetEnergyGlobalState(r, em);
  internal::gp.resize(gr.size());
  internal::gctd.resize(gr.size());
  internal::genergy = MAX_E;
  std::fill(internal::gp.begin(), internal::gp.end(), -1);
  std::fill(internal::gctd.begin(), internal::gctd.end(), CTD_NA);
  internal::gdp = array3d_t<energy_t, DP_SIZE>(gr.size() + 1);
  // TODO: Move this
  internal::gext = array2d_t<energy_t, EXT_SIZE>(gr.size() + 1);
}

}  // namespace mrna::mfe
