#include "common.h"
#include "parsing.h"
#include "constants.h"

namespace memerna {

namespace {
const int BUF_SIZE = 1024;
const uint32_t CRC_MAGIC = 0xEDB88320;

std::string SerialiseEnergyModel() {
  std::string data;

  // This isn't portable across machines with different endianness but I don't care.
#define APPEND_DATA(d) \
  do { \
    auto dp = reinterpret_cast<const char*>(&d); \
    data.insert(data.end(), dp, dp + sizeof(d)); \
  } while (0)

  APPEND_DATA(g_stack);
  APPEND_DATA(g_terminal);
  APPEND_DATA(g_internal_init);
  APPEND_DATA(g_internal_1x1);
  APPEND_DATA(g_internal_1x2);
  APPEND_DATA(g_internal_2x2);
  APPEND_DATA(g_internal_2x3_mismatch);
  APPEND_DATA(g_internal_other_mismatch);
  APPEND_DATA(g_internal_asym);
  APPEND_DATA(g_internal_augu_penalty);
  APPEND_DATA(g_internal_mismatch_1xk);
  APPEND_DATA(g_bulge_init);
  APPEND_DATA(g_bulge_special_c);
  APPEND_DATA(g_hairpin_init);
  APPEND_DATA(g_hairpin_uu_ga_first_mismatch);
  APPEND_DATA(g_hairpin_gg_first_mismatch);
  APPEND_DATA(g_hairpin_special_gu_closure);
  APPEND_DATA(g_hairpin_c3_loop);
  APPEND_DATA(g_hairpin_all_c_a);
  APPEND_DATA(g_hairpin_all_c_b);

  for (const auto& v : g_hairpin_e) {
    data += v.first;
    APPEND_DATA(v.second);
  }

  APPEND_DATA(g_multiloop_hack_a);
  APPEND_DATA(g_multiloop_hack_b);
  APPEND_DATA(g_dangle5_e);
  APPEND_DATA(g_dangle3_e);
  APPEND_DATA(g_coax_mismatch_non_contiguous);
  APPEND_DATA(g_coax_mismatch_wc_bonus);
  APPEND_DATA(g_coax_mismatch_gu_bonus);
  APPEND_DATA(g_augu_penalty);

  APPEND_DATA(constants::HAIRPIN_MIN_SZ);
  APPEND_DATA(constants::R);
  APPEND_DATA(constants::T);
  APPEND_DATA(constants::NINIO_MAX_ASYM);
  APPEND_DATA(constants::TWOLOOP_MAX_SZ);
#undef APPEND_DATA

  return data;
}

}

void LoadEnergyModelFromDataDir(const std::string& data_dir) {
  // Stacking interaction data.
  parsing::Parse2x2FromFile(data_dir + "/stacking.data", g_stack);

  // Terminal mismatch data.
  parsing::Parse2x2FromFile(data_dir + "/terminal.data", g_terminal);

  // Hairpin data.
  parsing::ParseMapFromFile(data_dir + "/hairpin.data", g_hairpin_e);
  parsing::ParseInitiationEnergyFromFile(data_dir + "/hairpin_initiation.data", g_hairpin_init);

  // Bulge loop data.
  parsing::ParseInitiationEnergyFromFile(data_dir + "/bulge_initiation.data", g_bulge_init);

  // Internal loop data.
  parsing::ParseInitiationEnergyFromFile(data_dir + "/internal_initiation.data", g_internal_init);
  parsing::ParseInternalLoop1x1FromFile(data_dir + "/internal_1x1.data");
  parsing::ParseInternalLoop1x2FromFile(data_dir + "/internal_1x2.data");
  parsing::ParseInternalLoop2x2FromFile(data_dir + "/internal_2x2.data");
  parsing::Parse2x2FromFile(data_dir + "/internal_2x3_mismatch.data", g_internal_2x3_mismatch);
  parsing::Parse2x2FromFile(data_dir + "/internal_other_mismatch.data", g_internal_other_mismatch);

  // Dangle data.
  parsing::ParseDangleDataFromFile(data_dir + "/dangle3.data", g_dangle3_e);
  parsing::ParseDangleDataFromFile(data_dir + "/dangle5.data", g_dangle5_e);

  // Other misc data.
  parsing::ParseMiscDataFromFile(data_dir + "/misc.data");
}

void LoadRandomEnergyModel(energy_t min_energy, energy_t max_energy) {
  verify_expr(min_energy <= max_energy, "Min energy must be <= max energy");
  static_assert(sizeof(energy_t) == 4, "Assumes size of energy_t is 4 bytes");
#define RANDOMISE_DATA(d) \
  do { \
    auto dp = reinterpret_cast<energy_t*>(&d); \
    for (unsigned int i = 0; i < sizeof(d) / sizeof(*dp); ++i) { \
      dp[i] = rand() % (max_energy - min_energy + 1) + min_energy; \
    } \
  } while (0)

  RANDOMISE_DATA(g_stack);
  RANDOMISE_DATA(g_terminal);
  RANDOMISE_DATA(g_internal_init);
  RANDOMISE_DATA(g_internal_1x1);
  RANDOMISE_DATA(g_internal_1x2);
  RANDOMISE_DATA(g_internal_2x2);
  RANDOMISE_DATA(g_internal_2x3_mismatch);
  RANDOMISE_DATA(g_internal_other_mismatch);
  RANDOMISE_DATA(g_internal_asym);
  RANDOMISE_DATA(g_internal_augu_penalty);
  RANDOMISE_DATA(g_internal_mismatch_1xk);
  RANDOMISE_DATA(g_bulge_init);
  RANDOMISE_DATA(g_bulge_special_c);
  RANDOMISE_DATA(g_hairpin_init);
  RANDOMISE_DATA(g_hairpin_uu_ga_first_mismatch);
  RANDOMISE_DATA(g_hairpin_gg_first_mismatch);
  RANDOMISE_DATA(g_hairpin_special_gu_closure);
  RANDOMISE_DATA(g_hairpin_c3_loop);
  RANDOMISE_DATA(g_hairpin_all_c_a);
  RANDOMISE_DATA(g_hairpin_all_c_b);

  for (auto& v : g_hairpin_e) {
    RANDOMISE_DATA(v.second);
  }

  RANDOMISE_DATA(g_multiloop_hack_a);
  RANDOMISE_DATA(g_multiloop_hack_b);
  RANDOMISE_DATA(g_dangle5_e);
  RANDOMISE_DATA(g_dangle3_e);
  RANDOMISE_DATA(g_coax_mismatch_non_contiguous);
  RANDOMISE_DATA(g_coax_mismatch_wc_bonus);
  RANDOMISE_DATA(g_coax_mismatch_gu_bonus);
  RANDOMISE_DATA(g_augu_penalty);
#undef RANDOMISE_DATA
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
  return Crc32(SerialiseEnergyModel());
}

}
