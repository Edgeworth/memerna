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
    folded_rna_t kFlushCoax = parsing::ParseViennaRna("GUGAAACACAAAAUGA", ".((...))((...)).");
    // NNDB T99 Multiloop example
    folded_rna_t kMultiloop = parsing::ParseViennaRna("UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");
};


TEST_F(EnergyTest, HairpinEnergy) {
  r = kHairpin1.r;
  EXPECT_EQ(AUGU_PENALTY + terminal_e[A][A][A][U] + HairpinInitiation(6), HairpinEnergy(3, 10));
  EXPECT_EQ(88, HairpinInitiation(87));
}

TEST_F(EnergyTest, BulgeLoopEnergy) {
  EXPECT_EQ(68, BulgeInitiation(57));
}

TEST_F(EnergyTest, InternalLoopEnergy) {
  EXPECT_EQ(46, InternalLoopInitiation(67));
}

TEST_F(EnergyTest, MultiloopEnergy) {
  EXPECT_EQ(57, MultiloopHackInitiation(4));
  EXPECT_EQ(74, MultiloopT99Initiation(8, 4));
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

TEST_F(EnergyTest, NNDBMultiloopExamples) {
  EXPECT_EQ(stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][A][U][G] +
            2 * AUGU_PENALTY + 2 * HairpinInitiation(3),
            ComputeEnergy(kFlushCoax));
  EXPECT_EQ(stacking_e[G][A][U][C] + terminal_e[C][G][A][G] + coax_mismatch_non_contiguous +
            3 * HairpinInitiation(3) + MultiloopT99Initiation(8, 4) + 2 * AUGU_PENALTY,
            ComputeEnergy(kMultiloop));
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
  EXPECT_EQ(
      AUGU_PENALTY * 2 + HairpinInitiation(3) +
      std::min(terminal_e[U][A][A][A], std::min(dangle3_e[U][A][A], dangle5_e[U][A][A])),
      ComputeEnergy(parsing::ParseViennaRna("AAAAAUA", ".(...).")));
  EXPECT_EQ(AUGU_PENALTY * 2 + HairpinInitiation(3), ComputeEnergy(parsing::ParseViennaRna("AAAAU", "(...)")));
}

}
}
