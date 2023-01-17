// Copyright 2016 Eliot Courtney.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "common_test.h"
#include "compute/energy/common/branch.h"
#include "compute/energy/common/t04like/branch.h"
#include "compute/energy/energy.h"
#include "compute/energy/t04/model.h"
#include "compute/energy/t04/precomp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"

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
    return test_t04_ems[GetParam()]
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
          em->terminal[A][A][A][U] + em->HairpinInitiation(6),
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][G][G][U] + em->hairpin_gg_first_mismatch + em->HairpinInitiation(5),
      GetEnergy(kNNDBHairpin2));

  if (em->hairpin.contains("CCGAGG")) {
    EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][C][G][G] +
            em->hairpin["CCGAGG"],
        GetEnergy(kNNDBHairpin3));
  }

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][C][C][U] + em->HairpinInitiation(6) + em->hairpin_all_c_a * 6 +
          em->hairpin_all_c_b,
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(em->stack[C][G][C][G] + em->stack[G][G][C][C] + em->stack[G][G][U][C] + em->gu_penalty +
          em->terminal[G][G][G][U] + em->hairpin_gg_first_mismatch + em->HairpinInitiation(5) +
          em->hairpin_special_gu_closure,
      GetEnergy(kNNDBHairpin5));
}

TEST_P(T22ModelTest, NNDBBulgeLoopExamples) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][C][G][G] + em->BulgeInitiation(1) +
          em->bulge_special_c + em->stack[C][G][C][G] + em->HairpinInitiation(3) -
          E(R * T * log(3)),
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(
      em->stack[G][A][U][C] + em->au_penalty + em->BulgeInitiation(3) + em->HairpinInitiation(3),
      GetEnergy(kNNDBBulge2));
}

TEST_P(T22ModelTest, NNDBInternalLoopExamples) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(5) +
          std::min(em->internal_asym, NINIO_MAX_ASYM) + em->internal_2x3_mismatch[A][G][G][U] +
          em->internal_2x3_mismatch[G][G][A][C] + em->internal_au_penalty +
          em->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] +
          em->internal_2x2[A][G][A][C][G][G][A][U] + em->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(6) +
          std::min(4 * em->internal_asym, NINIO_MAX_ASYM) + em->internal_au_penalty +
          em->HairpinInitiation(3),
      GetEnergy(kNNDBInternal1x5));
}

TEST_P(T22ModelTest, BaseCases) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->au_penalty + em->stack[G][A][U][C] + em->hairpin_init[3],
      GetEnergy("GAAAAUC", "((...))"));
  EXPECT_EQ(em->au_penalty + em->gu_penalty + em->stack[G][A][U][U] + em->hairpin_init[3],
      GetEnergy("GAAAAUU", "((...))"));
  EXPECT_EQ(em->au_penalty * 2 + em->HairpinInitiation(3) +
          std::min(ZERO_E,
              std::min(
                  em->terminal[U][A][A][A], std::min(em->dangle3[U][A][A], em->dangle5[U][A][A]))),
      GetEnergy("AAAAAUA", ".(...)."));
  EXPECT_EQ(em->au_penalty * 2 + em->HairpinInitiation(3), GetEnergy("AAAAU", "(...)"));
  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][U][A][G] + em->BulgeInitiation(1) +
          em->stack[U][G][C][A] + em->HairpinInitiation(3),
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

  // Special hairpin with helix affected by penultimate stack on one side only:
  EXPECT_EQ(E(1.51 + 0.22), GetEnergy(em, "AACAGUGCUU", "((......))"));
  EXPECT_EQ(E(0.81 + 0.22), GetEnergy(em, "AAACAGUGCUUU", ".((......))."));

  // Bulge loop:
  // GU/GU at 0, 1, 3, 4
  EXPECT_EQ(E(8.60 - 0.74 * 4), GetEnergy(em, "GGCGGAAAUUUU", "((.((...))))"));
  // GU/GU at 1, 2, 4, 5
  EXPECT_EQ(E(7.90 - 0.74 * 4), GetEnergy(em, "GGGCGGAAAUUUUU", ".((.((...))))."));

  // Bulge loop with lonely pairs:
  EXPECT_EQ(E(9.00), GetEnergy(em, "GCGAAAUU", "(.(...))"));
  EXPECT_EQ(E(8.30), GetEnergy(em, "GGCGAAAUUU", ".(.(...))."));

  // Internal loop:
  // AU/AU at 0, 1, 3, 4
  EXPECT_EQ(E(7.72 + 0.22 * 4), GetEnergy(em, "AAAAACCCUUCCCUU", "((.((...))...))"));
  // AU/AU at 1, 2, 4, 5
  EXPECT_EQ(E(7.02 + 0.22 * 4), GetEnergy(em, "AAAAAACCCUUCCCUUU", ".((.((...))...))."));

  // Internal loop with lonely pairs:
  EXPECT_EQ(E(9.60), GetEnergy(em, "AAACCCUCCCU", "(.(...)...)"));
  EXPECT_EQ(E(8.90), GetEnergy(em, "AAAACCCUCCCUU", ".(.(...)...)."));

  // Special internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU/AU at 0, 4
  EXPECT_EQ(E(6.92 + 0.22 * 2), GetEnergy(em, "AAAAACCCUUAUU", "((.((...)).))"));
  // AU/AU at 1, 5
  EXPECT_EQ(E(6.22 + 0.22 * 2), GetEnergy(em, "AAAAAACCCUUAUUU", ".((.((...)).))."));

  // Special internal loop with lonely pairs:
  EXPECT_EQ(E(8.80), GetEnergy(em, "AAACCCUAU", "(.(...).)"));
  EXPECT_EQ(E(8.10), GetEnergy(em, "AAAACCCUAUU", ".(.(...).)."));

  // Multiloop:
  EXPECT_EQ(E(23.70), GetEnergy(em, "ACACCCUCCACCCUCCACCCUCU", "(.(...)..(...)..(...).)"));
  EXPECT_EQ(E(23.00), GetEnergy(em, "AACACCCUCCACCCUCCACCCUCUU", ".(.(...)..(...)..(...).)."));

  // AU/AU at 0, 1
  EXPECT_EQ(
      E(22.76 + 0.22 * 2), GetEnergy(em, "AACACCCUCCACCCUCCACCCUCUU", "((.(...)..(...)..(...).))"));

  // AU/AU at 1, 2
  EXPECT_EQ(E(22.06 + 0.22 * 2),
      GetEnergy(em, "AAACACCCUCCACCCUCCACCCUCUUU", ".((.(...)..(...)..(...).))."));

  // AU/AU at 1, 2, 4, 5, 13, 14
  EXPECT_EQ(E(20.18 + 0.22 * 6),
      GetEnergy(em, "AAACAACCCUUCCAACCCUUCCACCCUCUUU", ".((.((...))..((...))..(...).))."));
}

#endif

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, T22ModelTest, testing::Range(0, NUM_TEST_MODELS));

}  // namespace mrna::erg
