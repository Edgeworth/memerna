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
const int RAND_MAX_HAIRPIN_SZ = 8;
const int RAND_MAX_NUM_HAIRPIN = 50;

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

void LoadRandomEnergyModel(uint32_t seed) {
  std::mt19937 eng(seed);
  std::uniform_int_distribution<memerna::energy_t> energy_dist(RAND_MIN_ENERGY, RAND_MAX_ENERGY);
  std::uniform_int_distribution<memerna::energy_t> nonneg_energy_dist(0, RAND_MAX_ENERGY);
#define RANDOMISE_DATA(d) \
  do { \
    auto dp = reinterpret_cast<energy_t*>(&d); \
    for (unsigned int i = 0; i < sizeof(d) / sizeof(*dp); ++i) { \
      dp[i] = energy_dist(eng); \
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
  g_internal_asym = nonneg_energy_dist(eng);  // This needs to be non-negative for some optimisations.
  RANDOMISE_DATA(g_internal_augu_penalty);
  RANDOMISE_DATA(g_bulge_init);
  RANDOMISE_DATA(g_bulge_special_c);
  RANDOMISE_DATA(g_hairpin_init);
  RANDOMISE_DATA(g_hairpin_uu_ga_first_mismatch);
  RANDOMISE_DATA(g_hairpin_gg_first_mismatch);
  RANDOMISE_DATA(g_hairpin_special_gu_closure);
  RANDOMISE_DATA(g_hairpin_c3_loop);
  RANDOMISE_DATA(g_hairpin_all_c_a);
  RANDOMISE_DATA(g_hairpin_all_c_b);

  g_hairpin_e.clear();
  std::uniform_int_distribution<int> hairpin_size_dist(constants::HAIRPIN_MIN_SZ, RAND_MAX_HAIRPIN_SZ);
  static_assert(constants::HAIRPIN_MIN_SZ <= RAND_MAX_HAIRPIN_SZ,
      "HAIRPIN_MIN_SZ > RAND_MAX_HAIRPIN does not make sense");
  std::uniform_int_distribution<int> num_hairpin_dist(constants::HAIRPIN_MIN_SZ, RAND_MAX_NUM_HAIRPIN);
  int num_hairpin = num_hairpin_dist(eng);
  for (int i = 0; i < num_hairpin; ++i) {
    auto hairpin = parsing::RnaToString(GenerateRandomRna(hairpin_size_dist(eng), eng));
    g_hairpin_e[hairpin] = energy_dist(eng);
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

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Correct things to be 180 degree rotations if required.
          g_stack[c][d][a][b] = g_stack[a][b][c][d];
          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              g_internal_1x1[d][e][f][a][b][c] = g_internal_1x1[a][b][c][d][e][f];
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  g_internal_2x2[e][f][g][h][a][b][c][d] = g_internal_2x2[a][b][c][d][e][f][g][h];
                }
              }
            }
          }
        }
      }
    }
  }

  std::string reason;
  verify_expr(IsValidEnergyModel(&reason), "invalid energy model: %s", reason.c_str());
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

  std::string reason;
  verify_expr(IsValidEnergyModel(&reason), "invalid energy model: %s", reason.c_str());
}

void LoadEnergyModelFromArgParse(const ArgParse& argparse) {
  if (argparse.HasFlag("seed")) {
    LoadRandomEnergyModel(uint32_t(atoi(argparse.GetOption("seed").c_str())));
  } else {
    LoadEnergyModelFromDataDir(argparse.GetOption("data-path"));
  }
}

bool IsValidEnergyModel(std::string* reason) {
#define CHECK_COND(cond, reason_str) \
  do { \
    if (!(cond)) { \
      if (reason) *reason = "expected " #cond ": " reason_str; \
      return false; \
    } \
  } while (0)

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Expect 180 degree rotations to be the same.
          CHECK_COND(g_stack[a][b][c][d] == g_stack[c][d][a][b],
              "180 degree rotations should be the same");
          CHECK_COND(g_internal_asym >= 0, "optimisations rely on this");

          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              CHECK_COND(g_internal_1x1[a][b][c][d][e][f] == g_internal_1x1[d][e][f][a][b][c],
                  "180 degree rotations should be the same");
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  CHECK_COND(g_internal_2x2[a][b][c][d][e][f][g][h] == g_internal_2x2[e][f][g][h][a][b][c][d],
                      "180 degree rotations should be the same");
                }
              }
            }
          }
        }
      }
    }
  }
#undef CHECK_COND
  return true;
}

}
}
