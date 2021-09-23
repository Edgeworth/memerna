// Copyright 2016 Eliot Courtney.
#include "common.h"

#include <random>

#include "parsing.h"

namespace mrna {

namespace {

const int BUF_SIZE = 1024;
const uint32_t CRC_MAGIC = 0xEDB88320;

}  // namespace

std::string sgetline(FILE* fp) {
  char buffer[BUF_SIZE];
  if (fgets(buffer, BUF_SIZE, fp) == nullptr) return "";
  const std::string s(buffer);
  // Minus one for null character.
  verify(s.size() < BUF_SIZE - 1, "buffer too small");
  return s;
}

std::string sfmt(const char* fmt, ...) {
  va_list l;
  va_start(l, fmt);
  std::string res = vsfmt(fmt, l);
  va_end(l);
  return res;
}

std::string vsfmt(const char* fmt, va_list l) {
  char buffer[BUF_SIZE];
  const int res = vsnprintf(buffer, BUF_SIZE, fmt, l);
  verify(res >= 0 && res < BUF_SIZE, "buffer too small");
  return buffer;
}

uint32_t Crc32(const std::string& data) {
  uint32_t table[1 << 8] = {};
  for (uint32_t i = 0; i < 1 << 8; ++i) {
    table[i] = i;
    for (int k = 0; k < 8; ++k)
      table[i] = (table[i] >> 1) ^ (static_cast<bool>(table[i] & 1) ? CRC_MAGIC : 0);
  }

  uint32_t window = 0xFFFFFFFF;
  for (char i : data) window = (window >> 8) ^ (table[(window & 0xFF) ^ uint8_t(i)]);

  return ~window;
}

}  // namespace mrna
