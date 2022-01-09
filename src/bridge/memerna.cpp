// Copyright 2016 Eliot Courtney.
#include "bridge/memerna.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/structure.h"

namespace mrna::bridge {

Energy Memerna::Efn(const Primary& r, const Secondary& s, std::string* desc) const {
  Computed computed;
  if (desc) {
    std::unique_ptr<energy::Structure> structure;
    computed = energy::ComputeEnergy(r, s, em_, &structure);
    for (const auto& s : structure->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    computed = energy::ComputeEnergy(r, s, em_);
  }

  return computed.energy;
}

Computed Memerna::Fold(const Primary& r) const {
  return std::get<Computed>(Context(r, em_, cfg_).Fold());
}

int Memerna::Suboptimal(
    subopt::SuboptimalCallback fn, const Primary& r, Energy energy_delta) const {
  return Context(r, em_, cfg_).Suboptimal(fn, true, energy_delta, -1);
}

std::vector<Computed> Memerna::SuboptimalIntoVector(const Primary& r, Energy energy_delta) const {
  return Context(r, em_, cfg_).SuboptimalIntoVector(true, energy_delta, -1);
}

std::pair<partition::Partition, Probabilities> Memerna::Partition(const Primary& r) const {
  auto partition = Context(r, em_, cfg_).Partition();
  auto probabilities = partition::ComputeProbabilities(partition);
  return {std::move(partition), std::move(probabilities)};
}

}  // namespace mrna::bridge
