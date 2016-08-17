#include <constants.h>
#include "fold/fold.h"
#include "energy/energy.h"
#include "parsing.h"
#include "gtest/gtest.h"
#include "common_test.h"

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
};


TEST_F(EnergyTest, MultiloopEnergy) {
  EXPECT_EQ(g_multiloop_hack_a + 4 * g_multiloop_hack_b, MultiloopInitiation(4));
}

TEST_F(EnergyTest, NNDBHairpinLoopExamples) {
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][A][U][G] + g_augu_penalty +
          g_terminal[A][A][A][U] + HairpinInitiation(6),
      ComputeEnergy(kNNDBHairpin1));
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][A][U][G] + g_augu_penalty +
          g_terminal[A][G][G][U] + g_hairpin_gg_first_mismatch + HairpinInitiation(5),
      ComputeEnergy(kNNDBHairpin2));
  EXPECT_EQ(g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][C][G][G] + g_hairpin_e["CCGAGG"],
      ComputeEnergy(kNNDBHairpin3));
  EXPECT_EQ(g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][A][U][G] +
      g_augu_penalty + g_terminal[A][C][C][U] + HairpinInitiation(6) + g_hairpin_all_c_a * 6 + g_hairpin_all_c_b,
      ComputeEnergy(kNNDBHairpin4));
  EXPECT_EQ(g_stack[C][G][C][G] + g_stack[G][G][C][C] + g_stack[G][G][U][C] + g_augu_penalty +
      g_terminal[G][G][G][U] + g_hairpin_gg_first_mismatch + HairpinInitiation(5) + g_hairpin_special_gu_closure,
      ComputeEnergy(kNNDBHairpin5));

  SetRna(kNNDBHairpin1.r);
  EXPECT_EQ(g_augu_penalty + g_terminal[A][A][A][U] + HairpinInitiation(6),
      fold::FastHairpin(3, 10, fold::PrecomputeFastHairpin()));
  SetRna(kNNDBHairpin2.r);
  EXPECT_EQ(
      g_augu_penalty + g_terminal[A][G][G][U] + g_hairpin_gg_first_mismatch + HairpinInitiation(5),
      fold::FastHairpin(3, 9, fold::PrecomputeFastHairpin()));
  SetRna(kNNDBHairpin3.r);
  EXPECT_EQ(g_hairpin_e["CCGAGG"], fold::FastHairpin(3, 8, fold::PrecomputeFastHairpin()));
  SetRna(kNNDBHairpin4.r);
  EXPECT_EQ(g_augu_penalty + g_terminal[A][C][C][U] + HairpinInitiation(6) +
      g_hairpin_all_c_a * 6 + g_hairpin_all_c_b,
      fold::FastHairpin(3, 10, fold::PrecomputeFastHairpin()));
  SetRna(kNNDBHairpin5.r);
  EXPECT_EQ(g_augu_penalty + g_terminal[G][G][G][U] + g_hairpin_gg_first_mismatch +
      HairpinInitiation(5) + g_hairpin_special_gu_closure,
      fold::FastHairpin(3, 9, fold::PrecomputeFastHairpin()));
}

TEST_F(EnergyTest, NNDBBulgeLoopExamples) {
  EXPECT_EQ(
      g_stack[G][C][G][C] + g_stack[C][C][G][G] + BulgeInitiation(1)
          + g_bulge_special_c + g_stack[C][G][C][G] + HairpinInitiation(3) -
          energy_t(round(10.0 * constants::R * constants::T * log(3))),
      ComputeEnergy(kNNDBBulge1));

  EXPECT_EQ(
      g_stack[G][A][U][C] + g_augu_penalty + BulgeInitiation(3) + HairpinInitiation(3),
      ComputeEnergy(kNNDBBulge2));
}


TEST_F(EnergyTest, NNDBMultiloopExamples) {
  EXPECT_EQ(g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][A][U][G] +
      2 * g_augu_penalty + 2 * HairpinInitiation(3),
      ComputeEnergy(kFlushCoax));
  EXPECT_EQ(g_stack[G][A][U][C] + g_terminal[C][G][A][G] + g_coax_mismatch_non_contiguous +
      3 * HairpinInitiation(3) + MultiloopInitiation(4) + 2 * g_augu_penalty,
      ComputeEnergy(kNNDBMultiloop));
}

TEST_F(EnergyTest, NNDBInternalLoopExamples) {
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[C][G][C][G] + InternalLoopInitiation(5) +
          std::min(g_internal_asym, constants::NINIO_MAX_ASYM) + g_internal_2x3_mismatch[A][G][G][U] +
          g_internal_2x3_mismatch[G][G][A][C] + g_internal_augu_penalty + HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x3));
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[C][G][C][G] + g_internal_2x2[A][G][A][C][G][G][A][U] +
          HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x2));
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[C][G][C][G] + InternalLoopInitiation(6) +
          std::min(4 * g_internal_asym, constants::NINIO_MAX_ASYM) + g_internal_augu_penalty + HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal1x5));
}

TEST_F(EnergyTest, BaseCases) {
  EXPECT_EQ(
      g_augu_penalty + g_stack[G][A][U][C] + g_hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketRna("GAAAAUC", "((...))")));
  EXPECT_EQ(
      g_augu_penalty * 2 + g_stack[G][A][U][U] + g_hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketRna("GAAAAUU", "((...))")));
  EXPECT_EQ(
      g_augu_penalty * 2 + HairpinInitiation(3) +
          std::min(g_terminal[U][A][A][A], std::min(g_dangle3_e[U][A][A], g_dangle5_e[U][A][A])),
      ComputeEnergy(parsing::ParseDotBracketRna("AAAAAUA", ".(...).")));
  EXPECT_EQ(g_augu_penalty * 2 + HairpinInitiation(3),
      ComputeEnergy(parsing::ParseDotBracketRna("AAAAU", "(...)")));
  EXPECT_EQ(
      g_stack[G][C][G][C] + g_stack[C][U][A][G] + BulgeInitiation(1) +
          g_stack[U][G][C][A] + HairpinInitiation(3),
      ComputeEnergy(kBulge1));
  EXPECT_EQ(
      InternalLoopInitiation(5) + g_internal_asym + g_internal_augu_penalty + g_augu_penalty +
          g_internal_2x3_mismatch[A][G][A][U] + g_internal_2x3_mismatch[C][A][A][G] + HairpinInitiation(3),
      ComputeEnergy(kInternal1));
}

TEST_F(EnergyTest, T04Tests) {
  if (EnergyModelChecksum() != T04_MODEL_HASH) return;

  EXPECT_EQ(88, HairpinInitiation(87));
  EXPECT_EQ(68, BulgeInitiation(57));
  EXPECT_EQ(46, InternalLoopInitiation(67));

  EXPECT_EQ(45, ComputeEnergy(parsing::ParseDotBracketRna("GCAAAGCC", "((...).)")));
  EXPECT_EQ(57, ComputeEnergy(parsing::ParseDotBracketRna("CCCAAAAUG", ".(.(...))")));
  EXPECT_EQ(55, ComputeEnergy(parsing::ParseDotBracketRna("UACAGA", "(....)")));
  EXPECT_EQ(-6, ComputeEnergy(parsing::ParseDotBracketRna("AGGGUCAUCCG", ".(((...))).")));
  EXPECT_EQ(80, ComputeEnergy(parsing::ParseDotBracketRna("AGAGAAACAAAU", "(..(...)...)")));
  EXPECT_EQ(95, ComputeEnergy(parsing::ParseDotBracketRna(
      "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)")));
  EXPECT_EQ(77, ComputeEnergy(parsing::ParseDotBracketRna("CCCGAAACAG", "(..(...).)")));
  EXPECT_EQ(74, ComputeEnergy(parsing::ParseDotBracketRna("GACAGAAACGCUGAAUC", "((..(...)......))")));
  EXPECT_EQ(173, ComputeEnergy(parsing::ParseDotBracketRna(
      "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))")));
  EXPECT_EQ(182, ComputeEnergy(parsing::ParseDotBracketRna(
      "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))")));
  EXPECT_EQ(176, ComputeEnergy(parsing::ParseDotBracketRna(
      "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)")));
  EXPECT_EQ(131, ComputeEnergy(parsing::ParseDotBracketRna(
      "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)")));
  EXPECT_EQ(-276, ComputeEnergy(parsing::ParseDotBracketRna(
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
      "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))....")));
  EXPECT_EQ(179, ComputeEnergy(parsing::ParseDotBracketRna(
      "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)")));

  // Special stacking - this is not implemented. TODO: Implement this?
  EXPECT_EQ(37, ComputeEnergy(parsing::ParseDotBracketRna("GGUCAAAGGUC", "((((...))))")));
  EXPECT_EQ(-45, ComputeEnergy(parsing::ParseDotBracketRna("GGGGAAACCCC", "((((...))))")));
  EXPECT_EQ(72, ComputeEnergy(parsing::ParseDotBracketRna("UGACAAAGGCGA", "(..(...)...)")));
}

}
}
