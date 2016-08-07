#ifndef MEMERNA_PARSING_H
#define MEMERNA_PARSING_H

#include <string>
#include <stack>
#include <unordered_map>
#include "base.h"
#include "globals.h"

namespace memerna {
namespace parsing {

rna_t ParseRnaFromString(const std::string& s);

folded_rna_t ParseDotBracketRna(const std::string& rna_str, const std::string& pairs_str);

std::string DotBracketFromPairs(const std::vector<int>& pairs);

std::string StringFromRna(const rna_t& rna);

void Parse2x2FromFile(const std::string& filename, energy_t (& output)[4][4][4][4]);

void ParseMapFromFile(const std::string& filename, std::unordered_map<std::string, energy_t>& output);

void ParseInitiationEnergyFromFile(const std::string& filename, energy_t (& output)[INITIATION_CACHE_SZ]);

void ParseInternalLoop1x1FromFile(const std::string& filename);

void ParseInternalLoop1x2FromFile(const std::string& filename);

void ParseInternalLoop2x2FromFile(const std::string& filename);

void ParseDangleDataFromFile(const std::string& filename, energy_t (& output)[4][4][4]);

void ParseMiscDataFromFile(const std::string& filename);

}
}

#endif //MEMERNA_PARSING_H
