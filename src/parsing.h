#ifndef MEMERNA_PARSING_H
#define MEMERNA_PARSING_H

#include <string>
#include <stack>
#include <unordered_map>
#include "base.h"
#include "globals.h"

namespace memerna {
namespace parsing {

primary_t StringToPrimary(const std::string& s);
std::string PrimaryToString(const primary_t& primary);
secondary_t ParseDotBracketSecondary(const std::string& prim_str, const std::string& pairs_str);
std::vector<int> DotBracketToPairs(const std::string& pairs_str);
std::string PairsToDotBracket(const std::vector<int>& pairs);
std::vector<int> BasePairListToPairs(const std::vector<std::pair<int, int>>& base_pairs, std::size_t size);


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
