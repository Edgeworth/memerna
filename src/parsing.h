#ifndef MEMERNA_PARSING_H
#define MEMERNA_PARSING_H

#include <string>
#include <stack>
#include <unordered_map>
#include "base.h"
#include "energy.h"
#include "globals.h"

namespace memerna {
namespace parsing {

rna_t ParseRnaFromString(const std::string& s);

folded_rna_t ParseViennaRna(const std::string& rna_str, const std::string& pairs_str);

void Parse2x2FromFile(const std::string& filename, energy::energy_t (& output)[4][4][4][4]);

void ParseMapFromFile(const std::string& filename, std::unordered_map<std::string, energy::energy_t>& output);

void ParseInitiationEnergyFromFile(const std::string& filename, energy::energy_t (& output)[INITIATION_CACHE_SZ]);

void ParseHairpinMiscDataFromFile(const std::string& filename);

void ParseBulgeMiscDataFromFile(const std::string& filename);

void ParseInternalLoop1x1FromFile(const std::string& filename);

void ParseInternalLoop1x2FromFile(const std::string& filename);

void ParseInternalLoop2x2FromFile(const std::string& filename);

void ParseInternalLoopMiscDataFromFile(const std::string& filename);

void ParseMultiloopMiscDataFromFile(const std::string& filename);

void ParseDangleDataFromFile(const std::string& filename, energy::energy_t (& output)[4][4][4]);

void ParseCoaxialMiscDataFromFile(const std::string& filename);

}
}

#endif //MEMERNA_PARSING_H
