// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_LOAD_MODEL_H_
#define COMPUTE_ENERGY_LOAD_MODEL_H_

#include <map>
#include <string>

#include "compute/energy/model.h"

namespace mrna::energy {

inline const std::map<std::string, Opt> ENERGY_OPTS = {
    {"seed", Opt("seed for random energy model for memerna").Arg()},
    {"rnastructure-data", Opt("data path for RNAstructure").Arg()},
    {"memerna-data", Opt("data path for given energy model for memerna").Arg()}};

EnergyModelPtr LoadEnergyModelFromDataDir(const std::string& data_dir);
EnergyModelPtr LoadRandomEnergyModel(uint_fast32_t seed);
EnergyModelPtr LoadEnergyModelFromArgParse(const ArgParse& args);

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_LOAD_MODEL_H_
