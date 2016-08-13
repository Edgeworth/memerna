#include <constants.h>
#include "energy.h"
#include "globals.h"
#include "parsing.h"
#include "gtest/gtest.h"

namespace memerna {
namespace energy {

class EnergyTest : public testing::Test {
public:
  folded_rna_t kNNDBHairpin1 = parsing::ParseDotBracketRna("CACAAAAAAAUGUG", "((((......))))");
  folded_rna_t kNNDBHairpin2 = parsing::ParseDotBracketRna("CACAGGAAGUGUG", "((((.....))))");
  folded_rna_t kNNDBHairpin3 = parsing::ParseDotBracketRna("CACCCGAGGGUG", "((((....))))");
  folded_rna_t kNNDBHairpin4 = parsing::ParseDotBracketRna("CACACCCCCCUGUG", "((((......))))");
  folded_rna_t kNNDBHairpin5 = parsing::ParseDotBracketRna("CGGGGGAAGUCCG", "((((.....))))");
  folded_rna_t kNNDBBulge1 = parsing::ParseDotBracketRna("GCCCGAAACGGC", "(((.(...))))");
  folded_rna_t kNNDBBulge2 = parsing::ParseDotBracketRna("GAACAGAAACUC", "((...(...)))");
  folded_rna_t kNNDBInternal2x3 = parsing::ParseDotBracketRna("CAGACGAAACGGAGUG", "((..((...))...))");
  folded_rna_t kNNDBInternal1x5 = parsing::ParseDotBracketRna("CAGCGAAACGGAAAGUG", "((.((...)).....))");
  folded_rna_t kNNDBInternal2x2 = parsing::ParseDotBracketRna("CAGACGAAACGGAUG", "((..((...))..))");
  folded_rna_t kFlushCoax = parsing::ParseDotBracketRna("GUGAAACACAAAAUGA", ".((...))((...)).");
  // NNDB T99 Multiloop example
  folded_rna_t kNNDBMultiloop = parsing::ParseDotBracketRna("UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");

  folded_rna_t kBulge1 = parsing::ParseDotBracketRna("GCUCGAAACAGC", "(((.(...))))");
  folded_rna_t kInternal1 = parsing::ParseDotBracketRna("AGAGAAACAAAU", "(..(...)...)");
  // 180 degree rotation of |kInternal1|.
  folded_rna_t kInternal2 = parsing::ParseDotBracketRna("CAAAUAAAAGAG", "(...(...)..)");
  folded_rna_t kInternal3 = parsing::ParseDotBracketRna("GGGGGAAACGGGC", "(...(...)...)");
};


TEST_F(EnergyTest, HairpinEnergy) {
  SetFoldedRna(kNNDBHairpin1);
  EXPECT_EQ(augu_penalty + terminal_e[A][A][A][U] + HairpinInitiation(6), HairpinEnergy(3, 10));
  EXPECT_EQ(88, HairpinInitiation(87));
}

TEST_F(EnergyTest, BulgeLoopEnergy) {
  EXPECT_EQ(68, BulgeInitiation(57));
}

TEST_F(EnergyTest, InternalLoopEnergy) {
  EXPECT_EQ(46, InternalLoopInitiation(67));
}

TEST_F(EnergyTest, MultiloopEnergy) {
  EXPECT_EQ(multiloop_hack_a + 4 * multiloop_hack_b, MultiloopHackInitiation(4));
  EXPECT_EQ(74, MultiloopT99Initiation(8, 4));
}

TEST_F(EnergyTest, NNDBHairpinLoopExamples) {
  EXPECT_EQ(
      stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][A][U][G] + augu_penalty +
      terminal_e[A][A][A][U] + HairpinInitiation(6),
      ComputeEnergy(kNNDBHairpin1));
  EXPECT_EQ(
      stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][A][U][G] + augu_penalty +
      terminal_e[A][G][G][U] + hairpin_gg_first_mismatch + HairpinInitiation(5),
      ComputeEnergy(kNNDBHairpin2));
  EXPECT_EQ(stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][C][G][G] + hairpin_e["CCGAGG"],
            ComputeEnergy(kNNDBHairpin3));
  EXPECT_EQ(stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][A][U][G] +
            augu_penalty + terminal_e[A][C][C][U] + HairpinInitiation(6) + hairpin_all_c_a * 6 + hairpin_all_c_b,
            ComputeEnergy(kNNDBHairpin4));
  EXPECT_EQ(stacking_e[C][G][C][G] + stacking_e[G][G][C][C] + stacking_e[G][G][U][C] + augu_penalty +
            terminal_e[G][G][G][U] + hairpin_gg_first_mismatch + HairpinInitiation(5) + hairpin_special_gu_closure,
            ComputeEnergy(kNNDBHairpin5));
}

TEST_F(EnergyTest, NNDBBulgeLoopExamples) {
  EXPECT_EQ(
      stacking_e[G][C][G][C] + stacking_e[C][C][G][G] + BulgeInitiation(1)
      + bulge_special_c + stacking_e[C][G][C][G] + HairpinInitiation(3) -
      energy_t(round(10.0 * constants::R * constants::T * log(3))),
      ComputeEnergy(kNNDBBulge1));

  EXPECT_EQ(
      stacking_e[G][A][U][C] + augu_penalty + BulgeInitiation(3) + HairpinInitiation(3),
      ComputeEnergy(kNNDBBulge2));
}


TEST_F(EnergyTest, NNDBMultiloopExamples) {
  EXPECT_EQ(stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][A][U][G] +
            2 * augu_penalty + 2 * HairpinInitiation(3),
            ComputeEnergy(kFlushCoax));
  EXPECT_EQ(stacking_e[G][A][U][C] + terminal_e[C][G][A][G] + coax_mismatch_non_contiguous +
            3 * HairpinInitiation(3) + MultiloopInitiation(8, 4) + 2 * augu_penalty,
            ComputeEnergy(kNNDBMultiloop));
}

TEST_F(EnergyTest, NNDBInternalLoopExamples) {
  EXPECT_EQ(
      stacking_e[C][A][U][G] + stacking_e[C][G][C][G] + InternalLoopInitiation(5) +
      std::min(1, constants::NINIO_MAX_ASYM) * internal_asym + internal_2x3_mismatch[A][G][G][U] +
      internal_2x3_mismatch[G][G][A][C] + internal_augu_penalty + HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x3));
  EXPECT_EQ(
      stacking_e[C][A][U][G] + stacking_e[C][G][C][G] + internal_2x2[A][G][A][C][G][G][A][U] +
      HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x2));
  EXPECT_EQ(
      stacking_e[C][A][U][G] + stacking_e[C][G][C][G] + InternalLoopInitiation(6) +
      std::min(4, constants::NINIO_MAX_ASYM) * internal_asym + internal_augu_penalty + HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal1x5));
}

TEST_F(EnergyTest, BaseCases) {
  EXPECT_EQ(
      augu_penalty + stacking_e[G][A][U][C] + hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketRna("GAAAAUC", "((...))")));
  EXPECT_EQ(
      augu_penalty * 2 + stacking_e[G][A][U][U] + hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketRna("GAAAAUU", "((...))")));
  EXPECT_EQ(
      augu_penalty * 2 + HairpinInitiation(3) +
      std::min(terminal_e[U][A][A][A], std::min(dangle3_e[U][A][A], dangle5_e[U][A][A])),
      ComputeEnergy(parsing::ParseDotBracketRna("AAAAAUA", ".(...).")));
  EXPECT_EQ(augu_penalty * 2 + HairpinInitiation(3), ComputeEnergy(parsing::ParseDotBracketRna("AAAAU", "(...)")));
  EXPECT_EQ(
      stacking_e[G][C][G][C] + stacking_e[C][U][A][G] + BulgeInitiation(1) +
      stacking_e[U][G][C][A] + HairpinInitiation(3),
      ComputeEnergy(kBulge1));

  EXPECT_EQ(
      InternalLoopInitiation(5) + internal_asym + internal_augu_penalty + augu_penalty +
      internal_2x3_mismatch[A][G][A][U] + internal_2x3_mismatch[C][A][A][G] + HairpinInitiation(3),
      ComputeEnergy(kInternal1));
}

}
}
