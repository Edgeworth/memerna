// Copyright 2016 Eliot Courtney.
#ifndef MODEL_PARSING_H_
#define MODEL_PARSING_H_

#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include "model/base.h"
#include "model/structure.h"

namespace mrna {

Primary StringToPrimary(const std::string& s);
std::string PrimaryToString(const Primary& r);
Secondary ParseDotBracketSecondary(const std::string& prim_str, const std::string& pairs_str);
std::vector<int> DotBracketToPairs(const std::string& pairs_str);
std::string PairsToDotBracket(const std::vector<int>& pairs);
Computed ParseCtdComputed(const std::string& prim_str, const std::string& pairs_str);
std::string ComputedToCtdString(const Computed& computed);
bool IsCtdString(const std::string& pairs_str);

}  // namespace mrna

#endif  // MODEL_PARSING_H_
