#include "bridge/memerna.h"
#include "energy/structure.h"

namespace memerna {
namespace bridge {

energy_t Memerna::Efn(const secondary_t& secondary, std::string* desc) const {
  computed_t computed;
  if (desc) {
    std::unique_ptr<energy::Structure> structure;
    computed = energy::ComputeEnergy(secondary, *em, &structure);
    for (const auto& s : structure->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    computed = energy::ComputeEnergy(secondary, *em);
  }

  return computed.energy;
}

computed_t Memerna::Fold(const primary_t& r) const { return fold::Context(r, em, options).Fold(); }

std::vector<computed_t> Memerna::Suboptimal(const primary_t& r, energy_t energy_delta) const {
  return fold::Context(r, em, options).SuboptimalSorted(energy_delta, -1);
}
}
}
