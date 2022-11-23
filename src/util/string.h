// Copyright 2021 Eliot Courtney.
#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include <cstdarg>
#include <cstdint>
#include <sstream>
#include <string>
#include <cctype>

namespace mrna {

std::string sgetline(std::istream& is);

std::string sfmt(const char* fmt, ...);

std::string vsfmt(const char* fmt, va_list l);

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

inline std::string ToLower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(c));
  return s;
}

std::string TrimLeft(const std::string& s);
std::string TrimRight(const std::string& s);
std::string Trim(const std::string& s);

uint32_t Crc32(const std::string& data);

}  // namespace mrna

#endif  // UTIL_STRING_H_
