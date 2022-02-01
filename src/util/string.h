// Copyright 2021 Eliot Courtney.
#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include <cstdint>
#include <cstdio>
#include <optional>
#include <sstream>
#include <string>

namespace mrna {

// Behaves as stated in the following cases:
//  error: throws exception
//  already eof: returns none
//  eof during: returns non-newline terminated line
//  normal: returns newline terminated line
std::optional<std::string> sgetline(FILE* fp);

std::string sfmt(const char* fmt, ...);

std::string vsfmt(const char* fmt, va_list l);

template <typename T>
T Conv(std::string s) {
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

uint32_t Crc32(const std::string& data);

}  // namespace mrna

#endif  // UTIL_STRING_H_
