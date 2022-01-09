// Copyright 2016 Eliot Courtney.
#include "bridge/memerna.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/structure.h"

namespace mrna::bridge {

energy::EnergyResult Memerna::Efn(const Primary& r, const Secondary& s, std::string* desc) const {
  energy::EnergyResult res;
  if (desc) {
    std::unique_ptr<energy::Structure> structure;
    res = energy::ComputeEnergy(r, s, nullptr, em_, &structure);
    for (const auto& s : structure->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    res = energy::ComputeEnergy(r, s, nullptr, em_);
  }

  return res;
}

FoldResult Memerna::Fold(const Primary& r) const { return Context(r, em_, cfg_).Fold(); }

int Memerna::Suboptimal(subopt::SuboptCallback fn, const Primary& r, Energy energy_delta) const {
  return Context(r, em_, cfg_).Suboptimal(fn, true, energy_delta, -1);
}

std::vector<subopt::SuboptResult> Memerna::SuboptimalIntoVector(
    const Primary& r, Energy energy_delta) const {
  return Context(r, em_, cfg_).SuboptimalIntoVector(true, energy_delta, -1);
}

partition::PartitionResult Memerna::Partition(const Primary& r) const {
  return Context(r, em_, cfg_).Partition();
}

}  // namespace mrna::bridge
