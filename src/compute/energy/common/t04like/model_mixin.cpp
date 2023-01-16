// Copyright 2023 Eliot Courtney.
#include "compute/energy/common/t04like/model_mixin.h"

#include "compute/energy/common/model.h"
#include "compute/energy/common/parse.h"

namespace mrna::erg {

bool T04ModelMixin::IsValid(std::string* reason) const {
  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Expect 180 degree rotations to be the same.
          CHECK_COND(
              stack[a][b][c][d] == stack[c][d][a][b], "180 degree rotations should be the same");
          CHECK_COND(internal_asym >= ZERO_E, "optimisations rely on this");

          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              CHECK_COND(internal_1x1[a][b][c][d][e][f] == internal_1x1[d][e][f][a][b][c],
                  "180 degree rotations should be the same");
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  CHECK_COND(
                      internal_2x2[a][b][c][d][e][f][g][h] == internal_2x2[e][f][g][h][a][b][c][d],
                      "180 degree rotations should be the same");
                }
              }
            }
          }
        }
      }
    }
  }
  return true;
}

void T04ModelMixin::ParseMiscDataFromFile(const std::string& filename) {
  std::ifstream f(filename);
  verify(f, "could not open file");

#define READ_DATA(var)                              \
  do {                                              \
    while (1) {                                     \
      auto line = sgetline(f);                      \
      verify(f, "unexpected EOF or error");         \
      if (line.empty() || line[0] == '/') continue; \
      (var) = Energy::FromString(Trim(line));       \
      break;                                        \
    }                                               \
  } while (0)

  // Bulge loops.
  READ_DATA(bulge_special_c);

  // Coaxial stacking.
  READ_DATA(coax_mismatch_non_contiguous);
  READ_DATA(coax_mismatch_wc_bonus);
  READ_DATA(coax_mismatch_gu_bonus);

  // Hairpin loops.
  READ_DATA(hairpin_uu_ga_first_mismatch);
  READ_DATA(hairpin_gg_first_mismatch);
  READ_DATA(hairpin_special_gu_closure);
  READ_DATA(hairpin_c3_loop);
  READ_DATA(hairpin_all_c_a);
  READ_DATA(hairpin_all_c_b);

  // Internal loops.
  READ_DATA(internal_asym);
  READ_DATA(internal_au_penalty);
  READ_DATA(internal_gu_penalty);

  // Multiloop data.
  READ_DATA(multiloop_hack_a);
  READ_DATA(multiloop_hack_b);

  // AU/GU penalty
  READ_DATA(au_penalty);
  READ_DATA(gu_penalty);
#undef READ_DATA

  f >> std::ws;
  verify(f.eof(), "expected EOF");
}

void T04ModelMixin::LoadFromDir(const std::string& data_dir) {
  // Stacking interaction data.
  Parse4MapFromFile(data_dir + "/stacking.data", stack);

  // Terminal mismatch data.
  Parse4MapFromFile(data_dir + "/terminal.data", terminal);

  // Hairpin data.
  ParseNMapFromFile(data_dir + "/hairpin.data", &hairpin);
  ParseVecFromFile(data_dir + "/hairpin_initiation.data", hairpin_init);

  // Bulge loop data.
  ParseVecFromFile(data_dir + "/bulge_initiation.data", bulge_init);

  // Internal loop data.
  ParseVecFromFile(data_dir + "/internal_initiation.data", internal_init);
  Parse6MapFromFile(data_dir + "/internal_1x1.data", internal_1x1);
  Parse7MapFromFile(data_dir + "/internal_1x2.data", internal_1x2);
  Parse8MapFromFile(data_dir + "/internal_2x2.data", internal_2x2);
  Parse4MapFromFile(data_dir + "/internal_2x3_mismatch.data", internal_2x3_mismatch);
  Parse4MapFromFile(data_dir + "/internal_other_mismatch.data", internal_other_mismatch);

  // Dangle data.
  Parse3MapFromFile(data_dir + "/dangle3.data", dangle3);
  Parse3MapFromFile(data_dir + "/dangle5.data", dangle5);

  // Other misc data.
  ParseMiscDataFromFile(data_dir + "/misc.data");
}

void T04ModelMixin::LoadRandom(std::mt19937& eng) {
  std::uniform_real_distribution<double> energy_dist(RAND_MIN_ENERGY, RAND_MAX_ENERGY);
  std::uniform_real_distribution<double> nonneg_energy_dist(0, RAND_MAX_ENERGY);

  RANDOMISE_DATA(stack);
  RANDOMISE_DATA(terminal);
  RANDOMISE_DATA(internal_init);
  RANDOMISE_DATA(internal_1x1);
  RANDOMISE_DATA(internal_1x2);
  RANDOMISE_DATA(internal_2x2);
  RANDOMISE_DATA(internal_2x3_mismatch);
  RANDOMISE_DATA(internal_other_mismatch);
  // This needs to be non-negative for some optimisations.
  internal_asym = Energy::FromDouble(nonneg_energy_dist(eng));
  RANDOMISE_DATA(internal_au_penalty);
  RANDOMISE_DATA(internal_gu_penalty);
  RANDOMISE_DATA(bulge_init);
  RANDOMISE_DATA(bulge_special_c);
  RANDOMISE_DATA(hairpin_init);
  RANDOMISE_DATA(hairpin_uu_ga_first_mismatch);
  RANDOMISE_DATA(hairpin_gg_first_mismatch);
  RANDOMISE_DATA(hairpin_special_gu_closure);
  RANDOMISE_DATA(hairpin_c3_loop);
  RANDOMISE_DATA(hairpin_all_c_a);
  RANDOMISE_DATA(hairpin_all_c_b);

  hairpin.clear();
  std::uniform_int_distribution<int> hairpin_size_dist(HAIRPIN_MIN_SZ, RAND_MAX_HAIRPIN_SZ);
  static_assert(HAIRPIN_MIN_SZ <= RAND_MAX_HAIRPIN_SZ,
      "HAIRPIN_MIN_SZ > RAND_MAX_HAIRPIN does not make sense");
  std::uniform_int_distribution<int> num_hairpin_dist(0, RAND_MAX_NUM_HAIRPIN);
  int num_hairpin = num_hairpin_dist(eng);
  for (int i = 0; i < num_hairpin; ++i) {
    auto hp = Primary::Random(hairpin_size_dist(eng)).ToSeq();
    hairpin[hp] = Energy::FromDouble(energy_dist(eng));
  }

  RANDOMISE_DATA(multiloop_hack_a);
  RANDOMISE_DATA(multiloop_hack_b);
  RANDOMISE_DATA(dangle5);
  RANDOMISE_DATA(dangle3);
  RANDOMISE_DATA(coax_mismatch_non_contiguous);
  RANDOMISE_DATA(coax_mismatch_wc_bonus);
  RANDOMISE_DATA(coax_mismatch_gu_bonus);
  RANDOMISE_DATA(au_penalty);
  RANDOMISE_DATA(gu_penalty);

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Correct things to be 180 degree rotations if required.
          stack[c][d][a][b] = stack[a][b][c][d];
          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              internal_1x1[d][e][f][a][b][c] = internal_1x1[a][b][c][d][e][f];
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  internal_2x2[e][f][g][h][a][b][c][d] = internal_2x2[a][b][c][d][e][f][g][h];
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace mrna::erg
