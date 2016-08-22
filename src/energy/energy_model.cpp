#include "energy/energy_model.h"
#include "parsing.h"
#include "common.h"
#include "constants.h"
#include <random>

namespace memerna {
namespace energy {

namespace {
const energy_t RAND_MIN_ENERGY = -1000;
const energy_t RAND_MAX_ENERGY = 1000;

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

void LoadRandomEnergyModel(int seed) {
  std::default_random_engine eng(seed);
  std::uniform_int_distribution<memerna::energy_t> uniform_dist(RAND_MIN_ENERGY, RAND_MAX_ENERGY);
  std::uniform_int_distribution<memerna::energy_t> uniform_dist_pos(0, RAND_MAX_ENERGY);
#define RANDOMISE_DATA(d) \
  do { \
    auto dp = reinterpret_cast<energy_t*>(&d); \
    for (unsigned int i = 0; i < sizeof(d) / sizeof(*dp); ++i) { \
      dp[i] = uniform_dist(eng); \
    } \
  } while (0);

  RANDOMISE_DATA(g_stack);
  RANDOMISE_DATA(g_terminal);
  RANDOMISE_DATA(g_internal_init);
  RANDOMISE_DATA(g_internal_1x1);
  RANDOMISE_DATA(g_internal_1x2);
  RANDOMISE_DATA(g_internal_2x2);
  RANDOMISE_DATA(g_internal_2x3_mismatch);
  RANDOMISE_DATA(g_internal_other_mismatch);
  g_internal_asym = uniform_dist_pos(eng);  // This needs to be non-negative for some optimisations.
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

uint32_t EnergyModelChecksum() {
  return Crc32(SerialiseEnergyModel());
}

void LoadEnergyModelFromDataDir(const std::string& data_dir) {
  // Stacking interaction data.
  memerna::parsing::Parse2x2FromFile(data_dir + "/stacking.data", memerna::g_stack);

  // Terminal mismatch data.
  memerna::parsing::Parse2x2FromFile(data_dir + "/terminal.data", memerna::g_terminal);

  // Hairpin data.
  memerna::parsing::ParseMapFromFile(data_dir + "/hairpin.data", memerna::g_hairpin_e);
  memerna::parsing::ParseInitiationEnergyFromFile(data_dir + "/hairpin_initiation.data", memerna::g_hairpin_init);

  // Bulge loop data.
  memerna::parsing::ParseInitiationEnergyFromFile(data_dir + "/bulge_initiation.data", memerna::g_bulge_init);

  // Internal loop data.
  memerna::parsing::ParseInitiationEnergyFromFile(data_dir + "/internal_initiation.data", memerna::g_internal_init);
  memerna::parsing::ParseInternalLoop1x1FromFile(data_dir + "/internal_1x1.data");
  memerna::parsing::ParseInternalLoop1x2FromFile(data_dir + "/internal_1x2.data");
  memerna::parsing::ParseInternalLoop2x2FromFile(data_dir + "/internal_2x2.data");
  memerna::parsing::Parse2x2FromFile(data_dir + "/internal_2x3_mismatch.data", memerna::g_internal_2x3_mismatch);
  memerna::parsing::Parse2x2FromFile(data_dir + "/internal_other_mismatch.data", memerna::g_internal_other_mismatch);

  // Dangle data.
  memerna::parsing::ParseDangleDataFromFile(data_dir + "/dangle3.data", memerna::g_dangle3_e);
  memerna::parsing::ParseDangleDataFromFile(data_dir + "/dangle5.data", memerna::g_dangle5_e);

  // Other misc data.
  memerna::parsing::ParseMiscDataFromFile(data_dir + "/misc.data");
}

void LoadEnergyModelFromArgParse(const ArgParse& argparse) {
  if (argparse.HasFlag("seed")) {
    LoadRandomEnergyModel(atoi(argparse.GetOption("seed").c_str()));
  } else {
    LoadEnergyModelFromDataDir(argparse.GetOption("data-path"));
  }
}

}
}
