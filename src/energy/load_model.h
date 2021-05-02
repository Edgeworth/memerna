// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#ifndef MEMERNA_LOAD_MODEL_H
#define MEMERNA_LOAD_MODEL_H

#include "common.h"
#include "energy/energy_model.h"

namespace memerna {
namespace energy {

const std::map<std::string, opt_t> ENERGY_OPTIONS = {
    {"seed", opt_t("seed for random energy model for memerna").Arg()},
    {"data-path", opt_t("data path for given energy model for memerna").Arg("data/")}};

EnergyModelPtr LoadEnergyModelFromDataDir(const std::string& data_dir);
EnergyModelPtr LoadRandomEnergyModel(uint_fast32_t seed);
EnergyModelPtr LoadEnergyModelFromArgParse(const ArgParse& argparse);
}  // namespace energy
}  // namespace memerna

#endif  // MEMERNA_LOAD_MODEL_H
