#ifndef MEMERNA_ENERGY_MODEL_H
#define MEMERNA_ENERGY_MODEL_H

#include <cstdarg>
#include <cstring>
#include "common.h"
#include "argparse.h"

namespace memerna {
namespace energy {

const std::map<std::string, ArgParse::option_t> ENERGY_OPTIONS = {
    {"seed", ArgParse::option_t("seed for random energy model for memerna").Arg()},
    {"data-path", ArgParse::option_t("data path for given energy model for memerna").Arg("data/")}
};

void LoadEnergyModelFromDataDir(const std::string& data_dir);
void LoadRandomEnergyModel(uint32_t seed);
void LoadEnergyModelFromArgParse(const ArgParse& argparse);
bool IsValidEnergyModel(std::string* reason = nullptr);
uint32_t EnergyModelChecksum();

}
}

#endif //MEMERNA_ENERGY_MODEL_H
