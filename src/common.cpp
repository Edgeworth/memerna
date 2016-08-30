#include <random>
#include "common.h"
#include "constants.h"
#include "parsing.h"

namespace memerna {

namespace {
const int BUF_SIZE = 1024;
const uint32_t CRC_MAGIC = 0xEDB88320;
}

std::string sgetline(FILE* fp) {
  char buffer[BUF_SIZE];
  if (fgets(buffer, BUF_SIZE, fp) == nullptr)
    return "";
  std::string s(buffer);
  // Minus one for null character.
  verify_expr(s.size() < BUF_SIZE - 1, "buffer too small");
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
  int res = vsnprintf(buffer, BUF_SIZE, fmt, l);
  verify_expr(res >= 0 && res < BUF_SIZE, "buffer too small");
  return buffer;
}

uint32_t Crc32(const std::string& data) {
  uint32_t table[1 << 8] = {};
  for (uint32_t i = 0; i < 1 << 8; ++i) {
    table[i] = i;
    for (int k = 0; k < 8; ++k) {
      table[i] = (table[i] >> 1) ^ ((table[i] & 1) ? CRC_MAGIC : 0);
    }
  }

  uint32_t window = 0xFFFFFFFF;
  for (uint32_t i = 0; i < data.size(); ++i) {
    window = (window >> 8) ^ (table[(window & 0xFF) ^ uint8_t(data[i])]);
  }

  return ~window;
}

computed_t::computed_t(const primary_t& primary)
  : r(primary), p(primary.size(), -1),
    ctds(primary.size(), CTD_NA), energy(constants::MAX_E) {}

computed_t::computed_t(const primary_t& primary, const std::vector<int>& p_,
    const std::vector<Ctd>& ctds_, energy_t energy_)
    : r(primary), p(p_), ctds(ctds_), energy(energy_) {
  assert(r.size() == p.size() && r.size() == ctds.size());
}
}
