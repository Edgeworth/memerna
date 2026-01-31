// Copyright 2021 Eliot Courtney.
#ifndef UTIL_ENUM_H_
#define UTIL_ENUM_H_

#include <fmt/core.h>

#include <boost/describe.hpp>
#include <cctype>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "util/error.h"

#define STRINGIFY(x) #x

#define MAKE_ENUM(name, ...) BOOST_DEFINE_ENUM_CLASS(name, __VA_ARGS__)
#define MAKE_NESTED_ENUM(name, ...) \
  enum class name { __VA_ARGS__ };  \
  BOOST_DESCRIBE_NESTED_ENUM(name, __VA_ARGS__)

namespace mrna {

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

template <typename T>
std::istream& operator>>(std::istream& is, T& value)
  requires boost::describe::has_describe_enumerators<T>::value
{
  std::string s;
  is >> s;
  s = NormalizeEnumName(s);
  for (const auto& [k, v] : EnumItems<T>()) {
    if (s == v) {
      value = k;
      return is;
    }
  }
  fatal("invalid enum value {}", s);
}

template <typename T>
std::ostream& operator<<(std::ostream& str, const T& value)
  requires boost::describe::has_describe_enumerators<T>::value
{
  auto name = boost::describe::enum_to_string(value, nullptr);
  verify(name != nullptr, "invalid enum value");
  return str << NormalizeEnumName(name);
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
    auto name = boost::describe::enum_to_string(value, nullptr);
    return fmt::formatter<fmt::string_view, char>().format(
        name ? mrna::NormalizeEnumName(name) : "?", ctx);
  }
};

#endif  // UTIL_ENUM_H_
