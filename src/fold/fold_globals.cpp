// Copyright 2016 E.
#include "fold/fold_globals.h"

#include "energy/energy_globals.h"
#include "globals.h"

namespace memerna {
namespace fold {
namespace internal {

std::vector<int> gp;
std::vector<Ctd> gctd;
std::string grep;
energy_t genergy;
array3d_t<energy_t, DP_SIZE> gdp;
array2d_t<energy_t, EXT_SIZE> gext;

}  // namespace internal

void SetFoldGlobalState(const primary_t& r, const energy::EnergyModel& em) {
  energy::SetEnergyGlobalState(r, em);
  internal::gp.resize(gr.size());
  internal::gctd.resize(gr.size());
  internal::genergy = MAX_E;
  std::fill(internal::gp.begin(), internal::gp.end(), -1);
  std::fill(internal::gctd.begin(), internal::gctd.end(), CTD_NA);
  internal::gdp = array3d_t<energy_t, internal::DP_SIZE>(gr.size() + 1);
  internal::gext = array2d_t<energy_t, internal::EXT_SIZE>(gr.size() + 1);
}

}  // namespace fold
}  // namespace memerna
