// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_ENERGY_H_
#define COMPUTE_ENERGY_ENERGY_H_

#include <memory>
#include <variant>

#include "model/constants.h"
#include "model/ctd.h"

namespace mrna::energy {

namespace t04 {
class Model;
using ModelPtr = std::shared_ptr<Model>;
}  // namespace t04

using EnergyModelPtr = std::variant<t04::ModelPtr>;

class Structure;

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

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_ENERGY_H_
