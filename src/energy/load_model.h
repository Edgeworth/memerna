// Copyright 2016 Eliot Courtney.
#ifndef ENERGY_LOAD_MODEL_H_
#define ENERGY_LOAD_MODEL_H_

#include <map>
#include <string>

#include "common.h"
#include "energy/energy_model.h"

namespace mrna::energy {

const std::map<std::string, opt_t> ENERGY_OPTIONS = {
    {"seed", opt_t("seed for random energy model for memerna").Arg()},
    {"rnastructure-data", opt_t("data path for RNAstructure").Arg()},
    {"memerna-data", opt_t("data path for given energy model for memerna").Arg()}};

EnergyModelPtr LoadEnergyModelFromDataDir(const std::string& data_dir);
EnergyModelPtr LoadRandomEnergyModel(uint_fast32_t seed);
EnergyModelPtr LoadEnergyModelFromArgParse(const ArgParse& args);

}  // namespace mrna::energy

#endif  // ENERGY_LOAD_MODEL_H_
