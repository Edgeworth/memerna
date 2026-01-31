// Copyright 2021 Eliot Courtney.
#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include <boost/algorithm/string/join.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "util/enum.h"  // IWYU pragma: keep - required for ADL to work for enum conversion.
#include "util/error.h"

namespace mrna {

std::string sgetline(std::istream& is);

template <typename T>
T Conv(const std::string& s) {
  T t{};
  std::stringstream ss(s);
  ss >> t;
  verify(ss && ss.eof(), "failed to convert '{}' to type", s);
  return t;
}

template <typename T>
std::string Conv(const T& s) {
  std::stringstream ss;
  ss << s;
  return ss.str();
}

std::string TrimLeft(const std::string& s);
std::string TrimRight(const std::string& s);
std::string Trim(const std::string& s);

std::vector<std::string> Split(const std::string& s, const std::string& delimiters);

template <typename Seq, typename Sep>
auto Join(const Seq& seq, const Sep& sep) {
  return boost::algorithm::join(seq, sep);
}

}  // namespace mrna

#endif  // UTIL_STRING_H_
