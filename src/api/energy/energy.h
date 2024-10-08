// Copyright 2022 E.
#ifndef API_ENERGY_ENERGY_H_
#define API_ENERGY_ENERGY_H_

#include <memory>

#include "model/ctd.h"
#include "model/energy.h"

namespace mrna {
class Structure;
}  // namespace mrna

namespace mrna::erg {

struct EnergyResult {
  EnergyResult() = default;
  EnergyResult(Energy energy, Ctds ctd, std::unique_ptr<Structure> struc);
  ~EnergyResult();

  EnergyResult(EnergyResult&&) = default;
  EnergyResult& operator=(EnergyResult&&) = default;
  EnergyResult(const EnergyResult&) = delete;
  EnergyResult& operator=(const EnergyResult&) = delete;

  Energy energy = ZERO_E;
  Ctds ctd;  // May be empty if CTDs were not computed.
  std::unique_ptr<Structure> struc;  // May be nullptr if you didn't ask for a Structure.
};

}  // namespace mrna::erg

#endif  // API_ENERGY_ENERGY_H_
