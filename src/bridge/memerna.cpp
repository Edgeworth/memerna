// Copyright 2016 E.
#include "bridge/memerna.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/structure.h"
#include "model/primary.h"

namespace mrna::bridge {

energy::EnergyResult Memerna::Efn(Primary r, Secondary s, std::string* desc) const {
  energy::EnergyResult res;
  if (desc) {
    res = ctx_.em().TotalEnergy(r, s, nullptr, true);
    for (const auto& s : res.struc->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    res = ctx_.em().TotalEnergy(r, s, nullptr, false);
  }

  return res;
}

FoldResult Memerna::Fold(Primary r) const { return ctx_.Fold(std::move(r)); }

int Memerna::Suboptimal(subopt::SuboptCallback fn, Primary r, Energy energy_delta) const {
  return ctx_.Suboptimal(std::move(r), fn, true, energy_delta, -1);
}

std::vector<subopt::SuboptResult> Memerna::SuboptimalIntoVector(
    Primary r, Energy energy_delta) const {
  return ctx_.SuboptimalIntoVector(std::move(r), true, energy_delta, -1);
}

partition::PartitionResult Memerna::Partition(Primary r) const {
  return ctx_.Partition(std::move(r));
}

Memerna Memerna::FromArgParse(const ArgParse& args) { return Memerna(Ctx::FromArgParse(args)); }

}  // namespace mrna::bridge
