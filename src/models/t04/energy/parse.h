#ifndef MODELS_T04_ENERGY_PARSE_H_
#define MODELS_T04_ENERGY_PARSE_H_

#include "models/common/boltz.h"
#include "models/common/model.h"
#include "models/common/parse.h"

namespace mrna::md::t04 {

template <typename EM>
void T04ParseMiscDataFromFile(EM& em, const std::string& filename) {
  std::ifstream f(filename);
  verify(f, "could not open file");

#define READ_DATA(var)                              \
  do {                                              \
    while (1) {                                     \
      auto line = sgetline(f);                      \
      verify(f, "unexpected EOF or error");         \
      if (line.empty() || line[0] == '/') continue; \
      em.var = Energy::FromString(Trim(line));      \
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

template <typename EM>
void T04LoadFromModelPath(EM& em, const std::string& path) {
  // Stacking interaction data.
  Parse4MapFromFile(path + "/stacking.data", em.stack);

  // Terminal mismatch data.
  Parse4MapFromFile(path + "/terminal.data", em.terminal);

  // Hairpin data.
  ParseNMapFromFile(path + "/hairpin.data", &em.hairpin);
  ParseVecFromFile(path + "/hairpin_initiation.data", em.hairpin_init);

  // Bulge loop data.
  ParseVecFromFile(path + "/bulge_initiation.data", em.bulge_init);

  // Internal loop data.
  ParseVecFromFile(path + "/internal_initiation.data", em.internal_init);
  Parse6MapFromFile(path + "/internal_1x1.data", em.internal_1x1);
  Parse7MapFromFile(path + "/internal_1x2.data", em.internal_1x2);
  Parse8MapFromFile(path + "/internal_2x2.data", em.internal_2x2);
  Parse4MapFromFile(path + "/internal_2x3_mismatch.data", em.internal_2x3_mismatch);
  Parse4MapFromFile(path + "/internal_other_mismatch.data", em.internal_other_mismatch);

  // Dangle data.
  Parse3MapFromFile(path + "/dangle3.data", em.dangle3);
  Parse3MapFromFile(path + "/dangle5.data", em.dangle5);

  // Other misc data.
  T04ParseMiscDataFromFile(em, path + "/misc.data");
}

template <typename EM>
void T04LoadRandom(EM& em, std::mt19937& eng, double min_energy, double max_energy,
    int max_hairpin_sz, int max_num_hairpin) {
  std::uniform_real_distribution<double> energy_dist(min_energy, max_energy);
  std::uniform_real_distribution<double> nonneg_energy_dist(0, max_energy);

  RANDOMISE_DATA(em, stack);
  RANDOMISE_DATA(em, terminal);
  RANDOMISE_DATA(em, internal_init);
  RANDOMISE_DATA(em, internal_1x1);
  RANDOMISE_DATA(em, internal_1x2);
  RANDOMISE_DATA(em, internal_2x2);
  RANDOMISE_DATA(em, internal_2x3_mismatch);
  RANDOMISE_DATA(em, internal_other_mismatch);
  // This needs to be non-negative for some optimisations.
  em.internal_asym = Energy::FromDouble(nonneg_energy_dist(eng));
  RANDOMISE_DATA(em, internal_au_penalty);
  RANDOMISE_DATA(em, internal_gu_penalty);
  RANDOMISE_DATA(em, bulge_init);
  RANDOMISE_DATA(em, bulge_special_c);
  RANDOMISE_DATA(em, hairpin_init);
  RANDOMISE_DATA(em, hairpin_uu_ga_first_mismatch);
  RANDOMISE_DATA(em, hairpin_gg_first_mismatch);
  RANDOMISE_DATA(em, hairpin_special_gu_closure);
  RANDOMISE_DATA(em, hairpin_c3_loop);
  RANDOMISE_DATA(em, hairpin_all_c_a);
  RANDOMISE_DATA(em, hairpin_all_c_b);

  em.hairpin.clear();
  std::uniform_int_distribution<int> hairpin_size_dist(HAIRPIN_MIN_SZ, max_hairpin_sz);
  verify(HAIRPIN_MIN_SZ <= max_hairpin_sz, "HAIRPIN_MIN_SZ > RAND_MAX_HAIRPIN does not make sense");
  std::uniform_int_distribution<int> num_hairpin_dist(0, max_num_hairpin);
  const int num_hairpin = num_hairpin_dist(eng);
  for (int i = 0; i < num_hairpin; ++i) {
    auto hp = Primary::Random(hairpin_size_dist(eng), eng).ToSeq();
    em.hairpin[hp] = Energy::FromDouble(energy_dist(eng));
  }

  RANDOMISE_DATA(em, multiloop_hack_a);
  RANDOMISE_DATA(em, multiloop_hack_b);
  RANDOMISE_DATA(em, dangle5);
  RANDOMISE_DATA(em, dangle3);
  RANDOMISE_DATA(em, coax_mismatch_non_contiguous);
  RANDOMISE_DATA(em, coax_mismatch_wc_bonus);
  RANDOMISE_DATA(em, coax_mismatch_gu_bonus);
  RANDOMISE_DATA(em, au_penalty);
  RANDOMISE_DATA(em, gu_penalty);

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Correct things to be 180 degree rotations if required.
          em.stack[c][d][a][b] = em.stack[a][b][c][d];
          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              em.internal_1x1[d][e][f][a][b][c] = em.internal_1x1[a][b][c][d][e][f];
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  em.internal_2x2[e][f][g][h][a][b][c][d] = em.internal_2x2[a][b][c][d][e][f][g][h];
                }
              }
            }
          }
        }
      }
    }
  }
}

template <typename EM>
bool T04IsValid(EM& em, std::string* reason) {
  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Expect 180 degree rotations to be the same.
          CHECK_COND(em.stack[a][b][c][d] == em.stack[c][d][a][b],
              "180 degree rotations should be the same");
          CHECK_COND(em.internal_asym >= ZERO_E, "optimisations rely on this");

          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              CHECK_COND(em.internal_1x1[a][b][c][d][e][f] == em.internal_1x1[d][e][f][a][b][c],
                  "180 degree rotations should be the same");
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  CHECK_COND(em.internal_2x2[a][b][c][d][e][f][g][h] ==
                          em.internal_2x2[e][f][g][h][a][b][c][d],
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

template <typename BEM>
void T04LoadBoltz(BEM& bem) {
  FILL_BOLTZ(bem, stack);
  FILL_BOLTZ(bem, terminal);
  FILL_BOLTZ(bem, internal_init);
  FILL_BOLTZ(bem, internal_1x1);
  FILL_BOLTZ(bem, internal_1x2);
  FILL_BOLTZ(bem, internal_2x2);
  FILL_BOLTZ(bem, internal_2x3_mismatch);
  FILL_BOLTZ(bem, internal_other_mismatch);
  FILL_BOLTZ(bem, internal_asym);
  FILL_BOLTZ(bem, internal_au_penalty);
  FILL_BOLTZ(bem, internal_gu_penalty);
  FILL_BOLTZ(bem, bulge_init);
  FILL_BOLTZ(bem, bulge_special_c);
  FILL_BOLTZ(bem, hairpin_init);
  FILL_BOLTZ(bem, hairpin_uu_ga_first_mismatch);
  FILL_BOLTZ(bem, hairpin_gg_first_mismatch);
  FILL_BOLTZ(bem, hairpin_special_gu_closure);
  FILL_BOLTZ(bem, hairpin_c3_loop);
  FILL_BOLTZ(bem, hairpin_all_c_a);
  FILL_BOLTZ(bem, hairpin_all_c_b);
  FILL_BOLTZ(bem, multiloop_hack_a);
  FILL_BOLTZ(bem, multiloop_hack_b);
  FILL_BOLTZ(bem, dangle5);
  FILL_BOLTZ(bem, dangle3);
  FILL_BOLTZ(bem, coax_mismatch_non_contiguous);
  FILL_BOLTZ(bem, coax_mismatch_wc_bonus);
  FILL_BOLTZ(bem, coax_mismatch_gu_bonus);
  FILL_BOLTZ(bem, au_penalty);
  FILL_BOLTZ(bem, gu_penalty);

  for (const auto& kv : bem.em().hairpin) bem.hairpin[kv.first] = kv.second.Boltz();
}

}  // namespace mrna::md::t04

#endif  // MODELS_T04_ENERGY_PARSE_H_
