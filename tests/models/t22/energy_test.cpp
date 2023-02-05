// Copyright 2016 Eliot Courtney.
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>

#include "common_test.h"
#include "gtest/gtest.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/t22/energy/model.h"

namespace mrna::erg {

Energy GetEnergy(const t22::Model::Ptr& em, const std::string& r, const std::string& db) {
  return em->TotalEnergy(Primary::FromSeq(r), Secondary::FromDb(db), nullptr).energy;
}

class T22ModelTest : public testing::TestWithParam<int> {
 public:
  std::tuple<Primary, Secondary> kNNDBHairpin1 = ParseSeqDb("CACAAAAAAAUGUG", "((((......))))");
  std::tuple<Primary, Secondary> kNNDBHairpin2 = ParseSeqDb("CACAGGAAGUGUG", "((((.....))))");
  std::tuple<Primary, Secondary> kNNDBHairpin3 = ParseSeqDb("CACCCGAGGGUG", "((((....))))");
  std::tuple<Primary, Secondary> kNNDBHairpin4 = ParseSeqDb("CACACCCCCCUGUG", "((((......))))");
  std::tuple<Primary, Secondary> kNNDBHairpin5 = ParseSeqDb("CGGGGGAAGUCCG", "((((.....))))");
  std::tuple<Primary, Secondary> kNNDBBulge1 = ParseSeqDb("GCCCGAAACGGC", "(((.(...))))");
  std::tuple<Primary, Secondary> kNNDBBulge2 = ParseSeqDb("GAACAGAAACUC", "((...(...)))");
  std::tuple<Primary, Secondary> kNNDBInternal2x3 =
      ParseSeqDb("CAGACGAAACGGAGUG", "((..((...))...))");
  std::tuple<Primary, Secondary> kNNDBInternal1x5 =
      ParseSeqDb("CAGCGAAACGGAAAGUG", "((.((...)).....))");
  std::tuple<Primary, Secondary> kNNDBInternal2x2 =
      ParseSeqDb("CAGACGAAACGGAUG", "((..((...))..))");

  std::tuple<Primary, Secondary> kBulge1 = ParseSeqDb("GCUCGAAACAGC", "(((.(...))))");
  std::tuple<Primary, Secondary> kInternal1 = ParseSeqDb("AGAGAAACAAAU", "(..(...)...)");

  static Energy GetEnergy(const std::string& r, const std::string& db) {
    return GetEnergy({Primary::FromSeq(r), Secondary::FromDb(db)});
  }

  static Energy GetEnergy(const std::tuple<Primary, Secondary>& s) {
    return test_t22_ems[GetParam()]
        ->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr)
        .energy;
  }
};

TEST_P(T22ModelTest, MultiloopEnergy) {
  auto em = test_t22_ems[GetParam()];
  EXPECT_EQ(em->multiloop_hack_a + 4 * em->multiloop_hack_b, em->MultiloopInitiation(4));
}

TEST_P(T22ModelTest, NNDBHairpinLoopExamples) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][A][A][U] + em->HairpinInitiation(6) + em->penultimate_stack[C][A][U][G] +
          em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][G][G][U] + em->hairpin_gg_first_mismatch + em->HairpinInitiation(5) +
          em->penultimate_stack[C][A][U][G] + em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBHairpin2));

  if (em->hairpin.contains("CCGAGG")) {
    EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][C][G][G] +
            em->hairpin["CCGAGG"] + em->penultimate_stack[C][C][G][G] +
            em->penultimate_stack[U][G][C][A],
        GetEnergy(kNNDBHairpin3));
  }

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][C][C][U] + em->HairpinInitiation(6) + em->hairpin_all_c_a * 6 +
          em->hairpin_all_c_b + em->penultimate_stack[C][A][U][G] +
          em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(em->stack[C][G][C][G] + em->stack[G][G][C][C] + em->stack[G][G][U][C] + em->gu_penalty +
          em->terminal[G][G][G][U] + em->hairpin_gg_first_mismatch + em->HairpinInitiation(5) +
          em->hairpin_special_gu_closure + em->penultimate_stack[G][G][U][C] +
          em->penultimate_stack[C][G][C][G],
      GetEnergy(kNNDBHairpin5));
}

TEST_P(T22ModelTest, NNDBBulgeLoopExamples) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][C][G][G] + em->BulgeInitiation(1) +
          em->bulge_special_c + em->stack[C][G][C][G] + em->HairpinInitiation(3) -
          E(R * T * log(3)) + em->penultimate_stack[C][G][C][G] + em->penultimate_stack[G][C][G][C],
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(em->stack[G][A][U][C] + em->au_penalty + em->BulgeInitiation(3) +
          em->HairpinInitiation(3) + em->penultimate_stack[U][C][G][A] +
          em->penultimate_stack[G][A][U][C],
      GetEnergy(kNNDBBulge2));
}

TEST_P(T22ModelTest, NNDBInternalLoopExamples) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(5) +
          std::min(em->internal_asym, NINIO_MAX_ASYM) + em->internal_2x3_mismatch[A][G][G][U] +
          em->internal_2x3_mismatch[G][G][A][C] + em->internal_au_penalty +
          em->HairpinInitiation(3) + em->penultimate_stack[C][G][C][G] +
          em->penultimate_stack[C][G][C][G] + em->penultimate_stack[C][A][U][G] +
          em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] +
          em->internal_2x2[A][G][A][C][G][G][A][U] + em->HairpinInitiation(3) +
          em->penultimate_stack[C][G][C][G] + em->penultimate_stack[C][G][C][G] +
          em->penultimate_stack[C][A][U][G] + em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(6) +
          std::min(4 * em->internal_asym, NINIO_MAX_ASYM) + em->internal_au_penalty +
          em->HairpinInitiation(3) + em->penultimate_stack[C][G][C][G] +
          em->penultimate_stack[C][G][C][G] + em->penultimate_stack[C][A][U][G] +
          em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBInternal1x5));
}

TEST_P(T22ModelTest, BaseCases) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->au_penalty + em->stack[G][A][U][C] + em->hairpin_init[3] +
          em->penultimate_stack[G][A][U][C] + em->penultimate_stack[U][C][G][A],
      GetEnergy("GAAAAUC", "((...))"));
  EXPECT_EQ(em->au_penalty + em->gu_penalty + em->stack[G][A][U][U] + em->hairpin_init[3] +
          em->penultimate_stack[G][A][U][U] + em->penultimate_stack[U][U][G][A],
      GetEnergy("GAAAAUU", "((...))"));
  EXPECT_EQ(em->au_penalty * 2 + em->HairpinInitiation(3) +
          std::min(ZERO_E,
              std::min(
                  em->terminal[U][A][A][A], std::min(em->dangle3[U][A][A], em->dangle5[U][A][A]))),
      GetEnergy("AAAAAUA", ".(...)."));
  EXPECT_EQ(em->au_penalty * 2 + em->HairpinInitiation(3), GetEnergy("AAAAU", "(...)"));
  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][U][A][G] + em->BulgeInitiation(1) +
          em->stack[U][G][C][A] + em->HairpinInitiation(3) + em->penultimate_stack[U][G][C][A] +
          em->penultimate_stack[G][C][G][C],
      GetEnergy(kBulge1));
  EXPECT_EQ(em->InternalLoopInitiation(5) + std::min(em->internal_asym, NINIO_MAX_ASYM) +
          em->internal_au_penalty + em->au_penalty + em->internal_2x3_mismatch[A][G][A][U] +
          em->internal_2x3_mismatch[C][A][A][G] + em->HairpinInitiation(3),
      GetEnergy(kInternal1));
}

#if ENERGY_PRECISION == 2

TEST(T22P2ModelTest, T22Tests) {
  auto em = t22p2;

  EXPECT_EQ(E(8.85), em->HairpinInitiation(87));
  EXPECT_EQ(E(6.79), em->BulgeInitiation(57));
  EXPECT_EQ(E(4.57), em->InternalLoopInitiation(67));

  // Example from https://doi.org/10.1093/nar/gkac261
  EXPECT_EQ(E(-3.65), GetEnergy(em, "UGUCGAUACCCUGUCGAUA", "((((((((...))))))))"));
  EXPECT_EQ(E(-4.05), GetEnergy(em, "UAGGUCAGCCCCUGGUCUA", "((((((((...))))))))"));

  // Penultimate stacking tests:

  // Hairpin + helices tests:
  // Test penultimate stacking not applied lonely pairs:
  EXPECT_EQ(E(5.70), GetEnergy(em, "AGCCCCU", ".(...)."));

  // AU end on CG (0.44):
  EXPECT_EQ(E(4.89 + 0.44), GetEnergy(em, "AGCCCCU", "((...))"));
  EXPECT_EQ(E(4.19 + 0.44), GetEnergy(em, "AAGCCCCUU", ".((...))."));

  // Counting twice, two AU on AUs (0.22):
  EXPECT_EQ(E(5.96 + 0.22 * 2), GetEnergy(em, "AACCCUU", "((...))"));
  EXPECT_EQ(E(5.26 + 0.22 * 2), GetEnergy(em, "AAACCCUUU", ".((...))."));

  // Special hairpin with lonely pair unaffected:
  EXPECT_EQ(E(2.45), GetEnergy(em, "ACAGUGCU", "(......)"));
  EXPECT_EQ(E(1.75), GetEnergy(em, "AAACAGUGCUUU", "..(......).."));

  // Special hairpins:
  // GU end on AU at 0; AU end on GU at 1
  EXPECT_EQ(E(1.85 - 0.31 - 0.71), GetEnergy(em, "GACAGUGCUU", "((......))"));
  // GU end on AU at 1; AU end on GU at 2
  EXPECT_EQ(E(1.15 - 0.31 - 0.71), GetEnergy(em, "AGACAGUGCUUU", ".((......))."));

  // Single nuc bulge loop - treated as continuous.
  // AU end on GU at 0, 4
  EXPECT_EQ(E(8.42 - 0.71 * 2), GetEnergy(em, "AGCGAAAAUUUU", "((.((...))))"));
  // AU end on GU at 1, 5
  EXPECT_EQ(E(7.72 - 0.71 * 2), GetEnergy(em, "GAGCGAAAAUUUUU", ".((.((...))))."));

  // Single nuc bulge, continuous single on outer:
  // GU end on GU at 0; AU end on GU at 3
  EXPECT_EQ(E(8.40 - 0.74 - 0.71), GetEnergy(em, "GCGAAAAUUU", "(.((...)))"));
  // GU end on GU at 1; AU end on GU at 4
  EXPECT_EQ(E(7.70 - 0.74 - 0.71), GetEnergy(em, "GGCGAAAAUUUU", ".(.((...)))."));

  // Single nuc bulge, continuous single on inner:
  // AU end on GU at 0, 4
  EXPECT_EQ(E(8.62 - 0.71 * 2), GetEnergy(em, "AGCAAAAUUU", "((.(...)))"));
  // AU end on GU at 1, 5
  EXPECT_EQ(E(7.92 - 0.71 * 2), GetEnergy(em, "GAGCAAAAUUUU", ".((.(...)))."));

  // Single nuc bulge, continuous single on both:
  // GU end on AU at 0; AU end on GU at 2
  EXPECT_EQ(E(8.60 - 0.31 - 0.71), GetEnergy(em, "GCAAAAUU", "(.(...))"));
  // GU end on AU at 1; AU end on GU at 3
  EXPECT_EQ(E(7.90 - 0.31 - 0.71), GetEnergy(em, "GGCAAAAUUU", ".(.(...))."));

  // Multi nuc bulge loop:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  EXPECT_EQ(E(7.62 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGCCGAAAAUUUU", "((..((...))))"));
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  EXPECT_EQ(E(7.54 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GAGCCAGAAAUUUUU", ".((..((...))))."));

  // Multi nuc bulge loop with lonely pairs:
  EXPECT_EQ(E(8.20), GetEnergy(em, "GCCGAAAUU", "(..(...))"));
  EXPECT_EQ(E(7.50), GetEnergy(em, "GGCCGAAAUUU", ".(..(...))."));

  // Internal loop:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  EXPECT_EQ(E(9.02 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGGGACCCUUCCCUU", "((.((...))...))"));
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  EXPECT_EQ(E(8.32 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GAGGGACCCUUCCCUUU", ".((.((...))...))."));

  // Internal loop with lonely pairs:
  EXPECT_EQ(E(9.60), GetEnergy(em, "AAACCCUCCCU", "(.(...)...)"));
  EXPECT_EQ(E(8.90), GetEnergy(em, "AAAACCCUCCCUU", ".(.(...)...)."));

  // Special 1x1 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  EXPECT_EQ(E(7.32 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGGGACCCUUCUU", "((.((...)).))"));
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  EXPECT_EQ(E(6.62 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GAGGGACCCUUCUUU", ".((.((...)).))."));

  // Special 1x1 internal loop with lonely pairs:
  EXPECT_EQ(E(7.90), GetEnergy(em, "AAACCCUAU", "(.(...).)"));
  EXPECT_EQ(E(7.20), GetEnergy(em, "AAAACCCUAUU", ".(.(...).)."));

  // TODO(0): Test 1x2, 2x2.

  // Multiloop:
  EXPECT_EQ(E(23.70), GetEnergy(em, "ACACCCUCCACCCUCCACCCUCU", "(.(...)..(...)..(...).)"));
  EXPECT_EQ(E(23.00), GetEnergy(em, "AACACCCUCCACCCUCCACCCUCUU", ".(.(...)..(...)..(...).)."));

  // AU end on GU at 0; GU end on AU at 1
  EXPECT_EQ(E(23.72 - 0.71 - 0.31),
      GetEnergy(em, "AGCACCCUCCACCCUCCACCCUCUU", "((.(...)..(...)..(...).))"));

  // AU end on GU at 1; GU end on AU at 2
  EXPECT_EQ(E(23.02 - 0.71 - 0.31),
      GetEnergy(em, "AAGCACCCUCCACCCUCCACCCUCUUU", ".((.(...)..(...)..(...).))."));

  // AU end on GU at 1, 4, 13; GU end on AU at 2, 5, 14
  EXPECT_EQ(E(23.06 - 0.71 * 3 - 0.31 * 3),
      GetEnergy(em, "AAGCAGCCCUUCCAGCCCUUCCACCCUCUUU", ".((.((...))..((...))..(...).))."));

  // Flush coax stack
  // GU end on AU at 0, 7; AU end on GU at 1, 8 (not taken as continuous)
  EXPECT_EQ(E(12.22 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GACCCUUGACCCUU", "((...))((...))"));

  // GU end on GC at 1; GU end on AU at 7; AU end on GU at 8
  EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31), GetEnergy(em, "GGCCCUCGACCCUU", "((...))((...))"));

  // GU end on AU at 1, 8; AU end on GU at 2, 9
  EXPECT_EQ(E(12.20 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGACCCUUGACCCUUU", ".((...))((...))."));

  // GU end on GC at 2; GU end on AU at 8; AU end on GU at 9
  EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31), GetEnergy(em, "AGGCCCUCGACCCUUU", ".((...))((...))."));

  // Flush coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  EXPECT_EQ(E(13.40), GetEnergy(em, "AGCCCUACCCUU", ".(...)(...)."));

  // Mismatch mediated coax stack
  // GU end on AU at 1, 9; AU end on GU at 2, 10
  EXPECT_EQ(E(9.40 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGACCCUUAGACCCUUU", ".((...)).((...))."));

  // GU end on GC at 2; GU end on AU at 9; AU end on GU at 10
  EXPECT_EQ(E(7.80 + 0.13 - 0.71 - 0.31), GetEnergy(em, "AGGCCCUCAGACCCUUU", ".((...)).((...))."));

  // Mismatch mediated coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  EXPECT_EQ(E(10.60), GetEnergy(em, "AGCCCUAACCCUU", ".(...).(...)."));
}

#endif

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, T22ModelTest, testing::Range(0, NUM_TEST_MODELS));

}  // namespace mrna::erg
