// Copyright 2016 Eliot Courtney.
#include "bridge/memerna.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/structure.h"

namespace mrna::bridge {

energy::EnergyResult Memerna::Efn(Primary r, Secondary s, std::string* desc) const {
  energy::EnergyResult res;
  if (desc) {
    std::unique_ptr<energy::Structure> structure;
    res = em_.TotalEnergy(r, s, nullptr, &structure);
    for (const auto& s : structure->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    res = em_.TotalEnergy(r, s, nullptr);
  }

  return res;
}

FoldResult Memerna::Fold(Primary r) const { return Context(em_, cfg_).Fold(std::move(r)); }

int Memerna::Suboptimal(subopt::SuboptCallback fn, Primary r, Energy energy_delta) const {
  return Context(em_, cfg_).Suboptimal(std::move(r), fn, true, energy_delta, -1);
}

std::vector<subopt::SuboptResult> Memerna::SuboptimalIntoVector(
    Primary r, Energy energy_delta) const {
  return Context(em_, cfg_).SuboptimalIntoVector(std::move(r), true, energy_delta, -1);
}

partition::PartitionResult Memerna::Partition(Primary r) const {
  return Context(em_, cfg_).Partition(std::move(r));
}

}  // namespace mrna::bridge
