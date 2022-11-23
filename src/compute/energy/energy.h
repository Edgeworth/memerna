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
class BoltzModel;
using BoltzModelPtr = std::shared_ptr<BoltzModel>;
}  // namespace t04

using EnergyModelPtr = std::variant<t04::ModelPtr>;
using BoltzEnergyModelPtr = std::variant<t04::BoltzModelPtr>;

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

inline bool CanPair(const EnergyModelPtr& em, const Primary& r, int st, int en) {
  return std::visit([&](const auto& em) { return em->CanPair(r, st, en); }, em);
}

inline EnergyResult TotalEnergy(const EnergyModelPtr& em, const Primary& r, const Secondary& s,
    const Ctds* given_ctd, bool build_structure = false) {
  return std::visit(
      [&](const auto& em) { return em->TotalEnergy(r, s, given_ctd, build_structure); }, em);
}

inline EnergyResult SubstructureEnergy(const EnergyModelPtr& em, const Primary& r,
    const Secondary& s, const Ctds* given_ctd, int st, int en, bool build_structure = false) {
  return std::visit(
      [&](const auto& em) {
        return em->SubstructureEnergy(r, s, given_ctd, st, en, build_structure);
      },
      em);
}

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_ENERGY_H_
