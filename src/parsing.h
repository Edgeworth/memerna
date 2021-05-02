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
#ifndef MEMERNA_PARSING_H
#define MEMERNA_PARSING_H

#include <stack>
#include <string>
#include <unordered_map>

#include "base.h"
#include "energy/energy_model.h"

namespace memerna {
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
}  // namespace memerna

#endif  // MEMERNA_PARSING_H
