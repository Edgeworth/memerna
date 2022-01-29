// Copyright 2022 E.
#ifndef COMPUTE_ENERGY_ENERGY_H_
#define COMPUTE_ENERGY_ENERGY_H_

#include <memory>

#include "model/ctd.h"
#include "model/model.h"

namespace mrna::energy {

class Structure;

struct EnergyResult {
  EnergyResult() = default;
  EnergyResult(Energy energy, Ctds ctd, std::unique_ptr<Structure> struc);
  ~EnergyResult();

  EnergyResult(EnergyResult&&) = default;
  EnergyResult& operator=(EnergyResult&&) = default;
  EnergyResult(const EnergyResult&) = delete;
  EnergyResult& operator=(const EnergyResult&) = delete;

  Energy energy = 0;
  Ctds ctd;  // May be empty if CTDs were not computed.
  std::unique_ptr<Structure> struc;  // May be nullptr if you didn't ask for a Structure.
};

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_ENERGY_H_
