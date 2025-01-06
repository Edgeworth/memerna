// Copyright 2021 Eliot Courtney.
#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include <boost/algorithm/string/join.hpp>
#include <boost/describe.hpp>
#include <cctype>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "util/error.h"

#define MAKE_ENUM(name, ...) BOOST_DEFINE_ENUM_CLASS(name, __VA_ARGS__)
#define MAKE_NESTED_ENUM(name, ...) \
  enum class name { __VA_ARGS__ };  \
  BOOST_DESCRIBE_NESTED_ENUM(name, __VA_ARGS__)

namespace mrna {

std::string sgetline(std::istream& is);

constexpr std::string NormalizeEnumName(std::string name) {
  for (char& c : name) {
    if (c == '_') c = '-';
    c = static_cast<char>(std::tolower(c));
  }
  return name;
}

template <typename T>
constexpr std::vector<T> EnumValues() {
  std::vector<T> values;
  boost::mp11::mp_for_each<boost::describe::describe_enumerators<T>>(
      [&](auto v) { values.push_back(v.value); });
  return values;
}

template <typename T>
constexpr std::vector<std::string> EnumNames() {
  std::vector<std::string> values;
  boost::mp11::mp_for_each<boost::describe::describe_enumerators<T>>(
      [&](auto v) { values.push_back(NormalizeEnumName(v.name)); });
  return values;
}

template <typename T>
constexpr std::vector<std::pair<T, std::string>> EnumItems() {
  std::vector<std::pair<T, std::string>> values;
  boost::mp11::mp_for_each<boost::describe::describe_enumerators<T>>(
      [&](auto v) { values.push_back({v.value, NormalizeEnumName(v.name)}); });
  return values;
}

template <typename T>
constexpr int EnumCount() {
  return boost::mp11::mp_size<boost::describe::describe_enumerators<T>>::value;
}

using boost::describe::operators::operator==;
using boost::describe::operators::operator!=;

template <typename T,
    typename std::enable_if_t<boost::describe::has_describe_enumerators<T>::value, bool> = true>
std::istream& operator>>(std::istream& is, T& value) {
  std::string s;
  is >> s;
  s = mrna::NormalizeEnumName(s);
  for (const auto& [k, v] : mrna::EnumItems<T>()) {
    if (s == v) {
      value = k;
      return is;
    }
  }
  fatal("invalid enum value {}", s);
}

template <typename T,
    typename std::enable_if_t<boost::describe::has_describe_enumerators<T>::value, bool> = true>
std::ostream& operator<<(std::ostream& str, const T& value) {
  auto name = boost::describe::enum_to_string(value, nullptr);
  verify(name != nullptr, "invalid enum value");
  return str << mrna::NormalizeEnumName(name);
}

template <typename T>
T Conv(const std::string& s) {
  T t;
  std::stringstream ss(s);
  ss >> t;
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

template <class T>
struct fmt::formatter<T, char,
    std::enable_if_t<boost::describe::has_describe_enumerators<T>::value>> {
 public:
  constexpr auto parse(format_parse_context& ctx) {
    return fmt::formatter<fmt::string_view, char>().parse(ctx);
  }

  auto format(const T& value, format_context& ctx) const {
    return fmt::formatter<fmt::string_view, char>().format(mrna::Conv(value), ctx);
  }
};

#endif  // UTIL_STRING_H_
