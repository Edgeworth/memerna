// Copyright 2016 E.
#include "bridge/memerna.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/structure.h"

namespace mrna::bridge {

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

computed_t Memerna::Fold(const primary_t& r) const { return Context(r, em, cfg).Fold(); }

int Memerna::Suboptimal(
    subopt::SuboptimalCallback fn, const primary_t& r, energy_t energy_delta) const {
  return Context(r, em, cfg).Suboptimal(fn, true, energy_delta, -1);
}

std::vector<computed_t> Memerna::SuboptimalIntoVector(
    const primary_t& r, energy_t energy_delta) const {
  return Context(r, em, cfg).SuboptimalIntoVector(true, energy_delta, -1);
}

std::pair<partition::partition_t, partition::probabilities_t> Memerna::Partition(
    const primary_t& r) const {
  auto partition = Context(r, em, cfg).Partition();
  auto probabilities = partition::ComputeProbabilities(partition);
  return {std::move(partition), std::move(probabilities)};
}

}  // namespace mrna::bridge
