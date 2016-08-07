#include "common.h"
#include "parsing.h"

namespace memerna {

namespace {
const int BUF_SIZE = 1024;
const uint32_t CRC_MAGIC = 0xEDB88320;
}


void Init() {
  // Stacking interaction data.
  parsing::Parse2x2FromFile("data/stacking.data", stacking_e);

  // Terminal mismatch data.
  parsing::Parse2x2FromFile("data/terminal.data", terminal_e);

  // Hairpin data.
  parsing::ParseMapFromFile("data/hairpin.data", hairpin_e);
  parsing::ParseInitiationEnergyFromFile("data/hairpin_initiation.data", hairpin_init);

  // Bulge loop data.
  parsing::ParseInitiationEnergyFromFile("data/bulge_initiation.data", bulge_init);

  // Internal loop data.
  parsing::ParseInitiationEnergyFromFile("data/internal_initiation.data", internal_init);
  parsing::ParseInternalLoop1x1FromFile("data/internal_1x1.data");
  parsing::ParseInternalLoop1x2FromFile("data/internal_1x2.data");
  parsing::ParseInternalLoop2x2FromFile("data/internal_2x2.data");
  parsing::Parse2x2FromFile("data/internal_2x3_mismatch.data", internal_2x3_mismatch);
  parsing::Parse2x2FromFile("data/internal_other_mismatch.data", internal_other_mismatch);

  // Dangle data.
  parsing::ParseDangleDataFromFile("data/dangle3.data", dangle3_e);
  parsing::ParseDangleDataFromFile("data/dangle5.data", dangle5_e);

  // Other misc data.
  parsing::ParseMiscDataFromFile("data/misc.data");
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


uint32_t EnergyModelChecksum() {
  std::string data;

#define APPEND_DATA(d) \
  do { \
    const char* dp = reinterpret_cast<const char*>(&d); \
    data.insert(data.end(), dp, dp + sizeof(d)); \
  } while (0)

  APPEND_DATA(stacking_e);
  APPEND_DATA(terminal_e);
  APPEND_DATA(internal_init);
  APPEND_DATA(internal_1x1);
  APPEND_DATA(internal_1x2);
  APPEND_DATA(internal_2x2);
  APPEND_DATA(internal_2x3_mismatch);
  APPEND_DATA(internal_other_mismatch);
  APPEND_DATA(internal_asym);
  APPEND_DATA(internal_augu_penalty);
  APPEND_DATA(internal_mismatch_1xk);
  APPEND_DATA(bulge_init);
  APPEND_DATA(bulge_special_c);
  APPEND_DATA(hairpin_init);
  APPEND_DATA(hairpin_uu_ga_first_mismatch);
  APPEND_DATA(hairpin_gg_first_mismatch);
  APPEND_DATA(hairpin_special_gu_closure);
  APPEND_DATA(hairpin_c3_loop);
  APPEND_DATA(hairpin_all_c_a);
  APPEND_DATA(hairpin_all_c_b);

  for (const auto& v : hairpin_e) {
    data += v.first;
    APPEND_DATA(v.second);
  }

  APPEND_DATA(multiloop_hack_a);
  APPEND_DATA(multiloop_hack_b);
  APPEND_DATA(multiloop_t99_a);
  APPEND_DATA(multiloop_t99_b);
  APPEND_DATA(multiloop_t99_c);
  APPEND_DATA(dangle5_e);
  APPEND_DATA(dangle3_e);
  APPEND_DATA(coax_mismatch_non_contiguous);
  APPEND_DATA(coax_mismatch_wc_bonus);
  APPEND_DATA(coax_mismatch_gu_bonus);
  APPEND_DATA(augu_penalty);
#undef APPEND_DATA

  return Crc32(data);
}

}
