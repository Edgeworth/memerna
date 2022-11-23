// Copyright 2022 Eliot Courtney.
#include <utility>

#include "compute/energy/model.h"
#include "compute/energy/structure.h"

namespace mrna::energy {

EnergyResult::EnergyResult(Energy energy, Ctds ctd, std::unique_ptr<Structure> struc)
    : energy(energy), ctd(std::move(ctd)), struc(std::move(struc)) {}

EnergyResult::~EnergyResult() = default;

}  // namespace mrna::energy
