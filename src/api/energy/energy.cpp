// Copyright 2022 Eliot Courtney.
#include "api/energy/energy.h"

#include <memory>
#include <utility>

#include "model/structure.h"

namespace mrna::erg {

EnergyResult::EnergyResult(Energy energy, Ctds ctd, std::unique_ptr<Structure> struc)
    : energy(energy), ctd(std::move(ctd)), struc(std::move(struc)) {}

EnergyResult::~EnergyResult() = default;

}  // namespace mrna::erg
