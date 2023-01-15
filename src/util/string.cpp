// Copyright 2016 Eliot Courtney.
#include "util/string.h"

#include <fmt/core.h>

#include <cctype>
#include <cstdarg>
#include <cstdio>
#include <memory>
#include <stdexcept>

namespace mrna {

namespace {

constexpr uint32_t CRC_MAGIC = 0xEDB88320;

}  // namespace

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

uint32_t Crc32(const std::string& data) {
  uint32_t table[1 << 8] = {};
  for (uint32_t i = 0; i < 1 << 8; ++i) {
    table[i] = i;
    for (int k = 0; k < 8; ++k)
      table[i] = (table[i] >> 1) ^ (static_cast<bool>(table[i] & 1) ? CRC_MAGIC : 0);
  }

  uint32_t window = 0xFFFFFFFF;
  for (char i : data) window = (window >> 8) ^ (table[(window & 0xFFU) ^ uint8_t(i)]);

  return ~window;
}

}  // namespace mrna
