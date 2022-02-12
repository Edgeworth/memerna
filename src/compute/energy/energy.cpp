// Copyright 2022 Eliot Courtney.
#include "compute/energy/energy.h"

#include <utility>

#include "compute/energy/structure.h"

namespace mrna::energy {

EnergyResult::EnergyResult(Energy energy, Ctds ctd, std::unique_ptr<Structure> struc)
    : energy(energy), ctd(std::move(ctd)), struc(std::move(struc)) {}

EnergyResult::~EnergyResult() = default;

}  // namespace mrna::energy
