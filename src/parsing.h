// Copyright 2016 E.
#ifndef PARSING_H_
#define PARSING_H_

#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include "base.h"
#include "energy/energy_model.h"

namespace mrna {
namespace parsing {

primary_t StringToPrimary(const std::string& s);
std::string PrimaryToString(const primary_t& r);
secondary_t ParseDotBracketSecondary(const std::string& prim_str, const std::string& pairs_str);
std::vector<int> DotBracketToPairs(const std::string& pairs_str);
std::string PairsToDotBracket(const std::vector<int>& pairs);
computed_t ParseCtdComputed(const std::string& prim_str, const std::string& pairs_str);
std::string ComputedToCtdString(const computed_t& computed);
bool IsCtdString(const std::string& pairs_str);
}  // namespace parsing
}  // namespace mrna

#endif  // PARSING_H_
