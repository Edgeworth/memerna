// Copyright 2016 Eliot Courtney.
#include "util/string.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/detail/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/type_index/type_index_facade.hpp>
#include <cctype>
#include <string>
#include <vector>

namespace mrna {

std::string sgetline(std::istream& is) {
  std::string s;
  std::getline(is, s);
  return s;
}

std::string TrimLeft(const std::string& s) {
  auto iter = s.begin();
  while (iter != s.end() && isspace(*iter)) ++iter;
  return {iter, s.end()};
}

std::string TrimRight(const std::string& s) {
  auto iter = s.end();
  while (iter != s.begin() && isspace(*(iter - 1))) --iter;
  return {s.begin(), iter};
}

std::string Trim(const std::string& s) { return TrimLeft(TrimRight(s)); }

std::vector<std::string> Split(const std::string& s, const std::string& delimiters) {
  std::vector<std::string> result;
  boost::split(result, s, boost::is_any_of(delimiters));
  return result;
}

}  // namespace mrna
