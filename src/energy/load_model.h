#ifndef MEMERNA_LOAD_MODEL_H
#define MEMERNA_LOAD_MODEL_H

#include "common.h"
#include "energy/energy_model.h"

namespace memerna {
namespace energy {

const std::map<std::string, ArgParse::option_t> ENERGY_OPTIONS = {
    {"seed", ArgParse::option_t("seed for random energy model for memerna").Arg()},
    {"data-path", ArgParse::option_t("data path for given energy model for memerna").Arg("data/")}};

EnergyModelPtr LoadEnergyModelFromDataDir(const std::string& data_dir);
EnergyModelPtr LoadRandomEnergyModel(uint_fast32_t seed);
EnergyModelPtr LoadEnergyModelFromArgParse(const ArgParse& argparse);
}
}

#endif  // MEMERNA_LOAD_MODEL_H
