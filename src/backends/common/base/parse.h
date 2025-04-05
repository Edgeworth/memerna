// Copyright 2025 Eliot Courtney.
#ifndef BACKENDS_COMMON_BASE_PARSE_H_
#define BACKENDS_COMMON_BASE_PARSE_H_

#include "backends/common/boltz.h"
#include "backends/common/model_mixin.h"
#include "backends/common/parse.h"

namespace mrna::md::base {

template <typename M>
void ParseMiscDataFromFile(M& m, const std::string& filename) {
  std::ifstream f(filename);
  verify(f, "could not open file");

#define READ_DATA(var)                              \
  do {                                              \
    while (1) {                                     \
      auto line = sgetline(f);                      \
      verify(f, "unexpected EOF or error");         \
      if (line.empty() || line[0] == '/') continue; \
      m.var = Energy::FromString(Trim(line));       \
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

template <typename M>
void LoadFromModelPath(M& m, const std::string& path) {
  // Stacking interaction data.
  Parse4MapFromFile(path + "/stacking.data", m.stack);

  // Terminal mismatch data.
  Parse4MapFromFile(path + "/terminal.data", m.terminal);

  // Hairpin data.
  ParseNMapFromFile(path + "/hairpin.data", &m.hairpin);
  ParseVecFromFile(path + "/hairpin_initiation.data", m.hairpin_init);

  // Bulge loop data.
  ParseVecFromFile(path + "/bulge_initiation.data", m.bulge_init);

  // Internal loop data.
  ParseVecFromFile(path + "/internal_initiation.data", m.internal_init);
  Parse6MapFromFile(path + "/internal_1x1.data", m.internal_1x1);
  Parse7MapFromFile(path + "/internal_1x2.data", m.internal_1x2);
  Parse8MapFromFile(path + "/internal_2x2.data", m.internal_2x2);
  Parse4MapFromFile(path + "/internal_2x3_mismatch.data", m.internal_2x3_mismatch);
  Parse4MapFromFile(path + "/internal_other_mismatch.data", m.internal_other_mismatch);

  // Dangle data.
  Parse3MapFromFile(path + "/dangle3.data", m.dangle3);
  Parse3MapFromFile(path + "/dangle5.data", m.dangle5);

  // Other misc data.
  ParseMiscDataFromFile(m, path + "/misc.data");
}

template <typename M>
void LoadRandomModel(M& m, std::mt19937& eng, double min_energy, double max_energy,
    int max_hairpin_sz, int max_num_hairpin) {
  std::uniform_real_distribution<double> energy_dist(min_energy, max_energy);
  std::uniform_real_distribution<double> nonneg_energy_dist(0, max_energy);

  RANDOMISE_DATA(m, stack);
  RANDOMISE_DATA(m, terminal);
  RANDOMISE_DATA(m, internal_init);
  RANDOMISE_DATA(m, internal_1x1);
  RANDOMISE_DATA(m, internal_1x2);
  RANDOMISE_DATA(m, internal_2x2);
  RANDOMISE_DATA(m, internal_2x3_mismatch);
  RANDOMISE_DATA(m, internal_other_mismatch);
  // This needs to be non-negative for some optimisations.
  m.internal_asym = Energy::FromFlt(nonneg_energy_dist(eng));
  RANDOMISE_DATA(m, internal_au_penalty);
  RANDOMISE_DATA(m, internal_gu_penalty);
  RANDOMISE_DATA(m, bulge_init);
  RANDOMISE_DATA(m, bulge_special_c);
  RANDOMISE_DATA(m, hairpin_init);
  RANDOMISE_DATA(m, hairpin_uu_ga_first_mismatch);
  RANDOMISE_DATA(m, hairpin_gg_first_mismatch);
  RANDOMISE_DATA(m, hairpin_special_gu_closure);
  RANDOMISE_DATA(m, hairpin_c3_loop);
  RANDOMISE_DATA(m, hairpin_all_c_a);
  RANDOMISE_DATA(m, hairpin_all_c_b);

  m.hairpin.clear();
  std::uniform_int_distribution<int> hairpin_size_dist(HAIRPIN_MIN_SZ, max_hairpin_sz);
  verify(HAIRPIN_MIN_SZ <= max_hairpin_sz, "HAIRPIN_MIN_SZ > RAND_MAX_HAIRPIN does not make sense");
  std::uniform_int_distribution<int> num_hairpin_dist(0, max_num_hairpin);
  const int num_hairpin = num_hairpin_dist(eng);
  for (int i = 0; i < num_hairpin; ++i) {
    auto hp = Primary::Random(hairpin_size_dist(eng), eng).ToSeq();
    m.hairpin[hp] = Energy::FromFlt(energy_dist(eng));
  }

  RANDOMISE_DATA(m, multiloop_hack_a);
  RANDOMISE_DATA(m, multiloop_hack_b);
  RANDOMISE_DATA(m, dangle5);
  RANDOMISE_DATA(m, dangle3);
  RANDOMISE_DATA(m, coax_mismatch_non_contiguous);
  RANDOMISE_DATA(m, coax_mismatch_wc_bonus);
  RANDOMISE_DATA(m, coax_mismatch_gu_bonus);
  RANDOMISE_DATA(m, au_penalty);
  RANDOMISE_DATA(m, gu_penalty);

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Correct things to be 180 degree rotations if required.
          m.stack[c][d][a][b] = m.stack[a][b][c][d];
          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              m.internal_1x1[d][e][f][a][b][c] = m.internal_1x1[a][b][c][d][e][f];
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  m.internal_2x2[e][f][g][h][a][b][c][d] = m.internal_2x2[a][b][c][d][e][f][g][h];
                }
              }
            }
          }
        }
      }
    }
  }
}

template <typename M>
bool ModelIsValid(M& m, std::string* reason) {
  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Expect 180 degree rotations to be the same.
          CHECK_COND(m.stack[a][b][c][d] == m.stack[c][d][a][b],
              "180 degree rotations should be the same");
          CHECK_COND(m.internal_asym >= ZERO_E, "optimisations rely on this");

          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              CHECK_COND(m.internal_1x1[a][b][c][d][e][f] == m.internal_1x1[d][e][f][a][b][c],
                  "180 degree rotations should be the same");
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  CHECK_COND(m.internal_2x2[a][b][c][d][e][f][g][h] ==
                          m.internal_2x2[e][f][g][h][a][b][c][d],
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

template <typename BM>
void LoadBoltzModel(BM& bm) {
  FILL_BOLTZ(bm, stack);
  FILL_BOLTZ(bm, terminal);
  FILL_BOLTZ(bm, internal_init);
  FILL_BOLTZ(bm, internal_1x1);
  FILL_BOLTZ(bm, internal_1x2);
  FILL_BOLTZ(bm, internal_2x2);
  FILL_BOLTZ(bm, internal_2x3_mismatch);
  FILL_BOLTZ(bm, internal_other_mismatch);
  FILL_BOLTZ(bm, internal_asym);
  FILL_BOLTZ(bm, internal_au_penalty);
  FILL_BOLTZ(bm, internal_gu_penalty);
  FILL_BOLTZ(bm, bulge_init);
  FILL_BOLTZ(bm, bulge_special_c);
  FILL_BOLTZ(bm, hairpin_init);
  FILL_BOLTZ(bm, hairpin_uu_ga_first_mismatch);
  FILL_BOLTZ(bm, hairpin_gg_first_mismatch);
  FILL_BOLTZ(bm, hairpin_special_gu_closure);
  FILL_BOLTZ(bm, hairpin_c3_loop);
  FILL_BOLTZ(bm, hairpin_all_c_a);
  FILL_BOLTZ(bm, hairpin_all_c_b);
  FILL_BOLTZ(bm, multiloop_hack_a);
  FILL_BOLTZ(bm, multiloop_hack_b);
  FILL_BOLTZ(bm, dangle5);
  FILL_BOLTZ(bm, dangle3);
  FILL_BOLTZ(bm, coax_mismatch_non_contiguous);
  FILL_BOLTZ(bm, coax_mismatch_wc_bonus);
  FILL_BOLTZ(bm, coax_mismatch_gu_bonus);
  FILL_BOLTZ(bm, au_penalty);
  FILL_BOLTZ(bm, gu_penalty);

  for (const auto& kv : bm.m().hairpin) bm.hairpin[kv.first] = kv.second.Boltz();
}

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_PARSE_H_
