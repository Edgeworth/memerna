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
T convert(std::string s) {
  T t;
  std::stringstream ss(s);
  ss >> t;
  return t;
}

uint32_t Crc32(const std::string& data);

}  // namespace mrna

#endif  // UTIL_STRING_H_
