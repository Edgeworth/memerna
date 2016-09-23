#include "globals.h"

namespace memerna {
namespace fold {
namespace internal {

primary_t gr;
std::vector<int> gp;
std::vector<Ctd> gctd;
precomp_t gpc;
energy_t genergy;
energy::EnergyModel gem;
array3d_t<energy_t, DP_SIZE> gdp;
array2d_t<energy_t, EXT_SIZE> gext;

void SetGlobalState(const primary_t& r, const energy::EnergyModel& em) {
  gr = r;
  gp.resize(gr.size());
  gctd.resize(gr.size());
  genergy = MAX_E;
  gem = em;
  gpc = PrecomputeData(gr, gem);
  std::fill(gp.begin(), gp.end(), -1);
  std::fill(gctd.begin(), gctd.end(), CTD_NA);
  gdp = array3d_t<energy_t, DP_SIZE>(gr.size() + 1);
  gext = array2d_t<energy_t, EXT_SIZE>(gr.size() + 1);
}
}
}
}
