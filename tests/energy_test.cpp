#include <energy.h>
#include <globals.h>
#include <parsing.h>
#include "gtest/gtest.h"

namespace memerna {
namespace energy {

class EnergyTest : public testing::Test {
public:
    folded_rna_t kHairpin1 = parsing::ParseViennaRna("CACAAAAAAAUGUG", "((((......))))");
    folded_rna_t kHairpin2 = parsing::ParseViennaRna("CACAGGAAGUGUG", "((((.....))))");
    folded_rna_t kHairpin3 = parsing::ParseViennaRna("CACCCGAGGGUG", "((((....))))");
    folded_rna_t kHairpin4 = parsing::ParseViennaRna("CACACCCCCCUGUG", "((((......))))");
    folded_rna_t kHairpin5 = parsing::ParseViennaRna("CGGGGGAAGUCCG", "((((.....))))");
    folded_rna_t kInternal2x3 = parsing::ParseViennaRna("CAGACGAAACGGAGUG", "((..((...))...))");
};


TEST_F(EnergyTest, HairpinEnergy) {
  r = kHairpin1.r;
  EXPECT_EQ(AUGU_PENALTY + terminal_e[A][A][A][U] + HairpinInitiation(6), HairpinEnergy(3, 10));
}

TEST_F(EnergyTest, NNDBHairpinLoopExamples) {
  EXPECT_EQ(
      stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][A][U][G] + AUGU_PENALTY +
      terminal_e[A][A][A][U] + HairpinInitiation(6),
      ComputeEnergy(kHairpin1));
  EXPECT_EQ(
      stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][A][U][G] + AUGU_PENALTY +
      terminal_e[A][G][G][U] + hairpin_gg_first_mismatch + HairpinInitiation(5),
      ComputeEnergy(kHairpin2));
  EXPECT_EQ(stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][C][G][G] + hairpin_e["CCGAGG"],
            ComputeEnergy(kHairpin3));
  EXPECT_EQ(stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][A][U][G] +
            AUGU_PENALTY + terminal_e[A][C][C][U] + HairpinInitiation(6) + hairpin_all_c_a * 6 + hairpin_all_c_b,
            ComputeEnergy(kHairpin4));
  EXPECT_EQ(stacking_e[C][G][C][G] + stacking_e[G][G][C][C] + stacking_e[G][G][U][C] + AUGU_PENALTY +
            terminal_e[G][G][G][U] + hairpin_gg_first_mismatch + HairpinInitiation(5) + hairpin_special_gu_closure,
            ComputeEnergy(kHairpin5));
}

TEST_F(EnergyTest, NNDBInternalLoopExamples) {
  EXPECT_EQ(
      stacking_e[C][A][U][G] + stacking_e[C][G][C][G] + InternalLoopInitiation(5) + internal_asym +
      internal_2x3_mismatch[A][G][G][U] + internal_2x3_mismatch[G][G][A][C] + internal_augu_penalty +
      HairpinInitiation(3),
      ComputeEnergy(kInternal2x3));
}

TEST_F(EnergyTest, BaseCases) {
  EXPECT_EQ(
      AUGU_PENALTY + stacking_e[G][A][U][C] + hairpin_init[3],
      ComputeEnergy(parsing::ParseViennaRna("GAAAAUC", "((...))")));
  EXPECT_EQ(
      AUGU_PENALTY * 2 + stacking_e[G][A][U][U] + hairpin_init[3],
      ComputeEnergy(parsing::ParseViennaRna("GAAAAUU", "((...))")));
  EXPECT_EQ(AUGU_PENALTY * 2 + HairpinInitiation(3), ComputeEnergy(parsing::ParseViennaRna("AAAAAUA", ".(...).")));
  EXPECT_EQ(AUGU_PENALTY * 2 + HairpinInitiation(3), ComputeEnergy(parsing::ParseViennaRna("AAAAU", "(...)")));
}

}
}
