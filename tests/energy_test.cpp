#include "constants.h"
#include "fold/fold.h"
#include "fold/fold_internal.h"
#include "parsing.h"
#include "gtest/gtest.h"
#include "common_test.h"
#include "energy/energy_model.h"

namespace memerna {
namespace energy {

class EnergyTest : public testing::Test {
public:
  secondary_t kNNDBHairpin1 = parsing::ParseDotBracketSecondary("CACAAAAAAAUGUG", "((((......))))");
  secondary_t kNNDBHairpin2 = parsing::ParseDotBracketSecondary("CACAGGAAGUGUG", "((((.....))))");
  secondary_t kNNDBHairpin3 = parsing::ParseDotBracketSecondary("CACCCGAGGGUG", "((((....))))");
  secondary_t kNNDBHairpin4 = parsing::ParseDotBracketSecondary("CACACCCCCCUGUG", "((((......))))");
  secondary_t kNNDBHairpin5 = parsing::ParseDotBracketSecondary("CGGGGGAAGUCCG", "((((.....))))");
  secondary_t kNNDBBulge1 = parsing::ParseDotBracketSecondary("GCCCGAAACGGC", "(((.(...))))");
  secondary_t kNNDBBulge2 = parsing::ParseDotBracketSecondary("GAACAGAAACUC", "((...(...)))");
  secondary_t kNNDBInternal2x3 = parsing::ParseDotBracketSecondary("CAGACGAAACGGAGUG", "((..((...))...))");
  secondary_t kNNDBInternal1x5 = parsing::ParseDotBracketSecondary("CAGCGAAACGGAAAGUG", "((.((...)).....))");
  secondary_t kNNDBInternal2x2 = parsing::ParseDotBracketSecondary("CAGACGAAACGGAUG", "((..((...))..))");
  secondary_t kFlushCoax = parsing::ParseDotBracketSecondary("GUGAAACACAAAAUGA", ".((...))((...)).");
  // NNDB T99 Multiloop example
  secondary_t kNNDBMultiloop = parsing::ParseDotBracketSecondary(
      "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");

  secondary_t kBulge1 = parsing::ParseDotBracketSecondary("GCUCGAAACAGC", "(((.(...))))");
  secondary_t kInternal1 = parsing::ParseDotBracketSecondary("AGAGAAACAAAU", "(..(...)...)");
};


TEST_F(EnergyTest, MultiloopEnergy) {
  EXPECT_EQ(g_multiloop_hack_a + 4 * g_multiloop_hack_b, MultiloopInitiation(4));
}

TEST_F(EnergyTest, NNDBHairpinLoopExamples) {
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][A][U][G] + g_augu_penalty +
          g_terminal[A][A][A][U] + HairpinInitiation(6),
      ComputeEnergy(kNNDBHairpin1).energy);
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][A][U][G] + g_augu_penalty +
          g_terminal[A][G][G][U] + g_hairpin_gg_first_mismatch + HairpinInitiation(5),
      ComputeEnergy(kNNDBHairpin2).energy);
  EXPECT_EQ(g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][C][G][G] + g_hairpin["CCGAGG"],
      ComputeEnergy(kNNDBHairpin3).energy);
  EXPECT_EQ(g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][A][U][G] +
      g_augu_penalty + g_terminal[A][C][C][U] + HairpinInitiation(6) + g_hairpin_all_c_a * 6 + g_hairpin_all_c_b,
      ComputeEnergy(kNNDBHairpin4).energy);
  EXPECT_EQ(g_stack[C][G][C][G] + g_stack[G][G][C][C] + g_stack[G][G][U][C] + g_augu_penalty +
      g_terminal[G][G][G][U] + g_hairpin_gg_first_mismatch + HairpinInitiation(5) + g_hairpin_special_gu_closure,
      ComputeEnergy(kNNDBHairpin5).energy);

  using fold::internal::FastHairpin;
  using fold::internal::PrecomputeFastHairpin;
  EXPECT_EQ(g_augu_penalty + g_terminal[A][A][A][U] + HairpinInitiation(6),
      FastHairpin(kNNDBHairpin1.r, 3, 10, PrecomputeFastHairpin(kNNDBHairpin1.r)));
  EXPECT_EQ(
      g_augu_penalty + g_terminal[A][G][G][U] + g_hairpin_gg_first_mismatch + HairpinInitiation(5),
      FastHairpin(kNNDBHairpin2.r, 3, 9, PrecomputeFastHairpin(kNNDBHairpin2.r)));
  EXPECT_EQ(g_hairpin["CCGAGG"],
      FastHairpin(kNNDBHairpin3.r, 3, 8, PrecomputeFastHairpin(kNNDBHairpin3.r)));
  EXPECT_EQ(g_augu_penalty + g_terminal[A][C][C][U] + HairpinInitiation(6) +
      g_hairpin_all_c_a * 6 + g_hairpin_all_c_b,
      FastHairpin(kNNDBHairpin4.r, 3, 10, PrecomputeFastHairpin(kNNDBHairpin4.r)));
  EXPECT_EQ(g_augu_penalty + g_terminal[G][G][G][U] + g_hairpin_gg_first_mismatch +
      HairpinInitiation(5) + g_hairpin_special_gu_closure,
      FastHairpin(kNNDBHairpin5.r, 3, 9, PrecomputeFastHairpin(kNNDBHairpin5.r)));
}

TEST_F(EnergyTest, NNDBBulgeLoopExamples) {
  EXPECT_EQ(
      g_stack[G][C][G][C] + g_stack[C][C][G][G] + BulgeInitiation(1)
          + g_bulge_special_c + g_stack[C][G][C][G] + HairpinInitiation(3) -
          energy_t(round(10.0 * constants::R * constants::T * log(3))),
      ComputeEnergy(kNNDBBulge1).energy);

  EXPECT_EQ(
      g_stack[G][A][U][C] + g_augu_penalty + BulgeInitiation(3) + HairpinInitiation(3),
      ComputeEnergy(kNNDBBulge2).energy);
}


TEST_F(EnergyTest, NNDBMultiloopExamples) {
  EXPECT_EQ(g_stack[C][A][U][G] + g_stack[A][C][G][U] + g_stack[C][A][U][G] +
      2 * g_augu_penalty + 2 * HairpinInitiation(3),
      ComputeEnergy(kFlushCoax).energy);
  EXPECT_EQ(g_stack[G][A][U][C] + g_terminal[C][G][A][G] + g_coax_mismatch_non_contiguous +
      3 * HairpinInitiation(3) + MultiloopInitiation(4) + 2 * g_augu_penalty,
      ComputeEnergy(kNNDBMultiloop).energy);
}

TEST_F(EnergyTest, NNDBInternalLoopExamples) {
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[C][G][C][G] + InternalLoopInitiation(5) +
          std::min(g_internal_asym, constants::NINIO_MAX_ASYM) + g_internal_2x3_mismatch[A][G][G][U] +
          g_internal_2x3_mismatch[G][G][A][C] + g_internal_augu_penalty + HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x3).energy);
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[C][G][C][G] + g_internal_2x2[A][G][A][C][G][G][A][U] +
          HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x2).energy);
  EXPECT_EQ(
      g_stack[C][A][U][G] + g_stack[C][G][C][G] + InternalLoopInitiation(6) +
          std::min(4 * g_internal_asym, constants::NINIO_MAX_ASYM) + g_internal_augu_penalty + HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal1x5).energy);
}

TEST_F(EnergyTest, BaseCases) {
  EXPECT_EQ(
      g_augu_penalty + g_stack[G][A][U][C] + g_hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketSecondary("GAAAAUC", "((...))")).energy);
  EXPECT_EQ(
      g_augu_penalty * 2 + g_stack[G][A][U][U] + g_hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketSecondary("GAAAAUU", "((...))")).energy);
  EXPECT_EQ(
      g_augu_penalty * 2 + HairpinInitiation(3) +
          std::min(g_terminal[U][A][A][A], std::min(g_dangle3[U][A][A], g_dangle5[U][A][A])),
      ComputeEnergy(parsing::ParseDotBracketSecondary("AAAAAUA", ".(...).")).energy);
  EXPECT_EQ(g_augu_penalty * 2 + HairpinInitiation(3),
      ComputeEnergy(parsing::ParseDotBracketSecondary("AAAAU", "(...)")).energy);
  EXPECT_EQ(
      g_stack[G][C][G][C] + g_stack[C][U][A][G] + BulgeInitiation(1) +
          g_stack[U][G][C][A] + HairpinInitiation(3),
      ComputeEnergy(kBulge1).energy);
  EXPECT_EQ(
      InternalLoopInitiation(5) + g_internal_asym + g_internal_augu_penalty + g_augu_penalty +
          g_internal_2x3_mismatch[A][G][A][U] + g_internal_2x3_mismatch[C][A][A][G] + HairpinInitiation(3),
      ComputeEnergy(kInternal1).energy);
}

TEST_F(EnergyTest, T04Tests) {
  ONLY_FOR_THIS_MODEL(T04_MODEL_HASH);

  EXPECT_EQ(88, HairpinInitiation(87));
  EXPECT_EQ(68, BulgeInitiation(57));
  EXPECT_EQ(46, InternalLoopInitiation(67));

  EXPECT_EQ(45, ComputeEnergy(parsing::ParseDotBracketSecondary("GCAAAGCC", "((...).)")).energy);
  EXPECT_EQ(57, ComputeEnergy(parsing::ParseDotBracketSecondary("CCCAAAAUG", ".(.(...))")).energy);
  EXPECT_EQ(55, ComputeEnergy(parsing::ParseDotBracketSecondary("UACAGA", "(....)")).energy);
  EXPECT_EQ(-6, ComputeEnergy(parsing::ParseDotBracketSecondary("AGGGUCAUCCG", ".(((...))).")).energy);
  EXPECT_EQ(80, ComputeEnergy(parsing::ParseDotBracketSecondary("AGAGAAACAAAU", "(..(...)...)")).energy);
  EXPECT_EQ(95, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)")).energy);
  EXPECT_EQ(77, ComputeEnergy(parsing::ParseDotBracketSecondary("CCCGAAACAG", "(..(...).)")).energy);
  EXPECT_EQ(74, ComputeEnergy(parsing::ParseDotBracketSecondary("GACAGAAACGCUGAAUC", "((..(...)......))")).energy);
  EXPECT_EQ(173, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))")).energy);
  EXPECT_EQ(182, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))")).energy);
  EXPECT_EQ(176, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)")).energy);
  EXPECT_EQ(131, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)")).energy);
  EXPECT_EQ(-276, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
      "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))....")).energy);
  EXPECT_EQ(179, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)")).energy);

  // Special stacking - this is not implemented. TODO: Implement this?
  EXPECT_EQ(37, ComputeEnergy(parsing::ParseDotBracketSecondary("GGUCAAAGGUC", "((((...))))")).energy);
  EXPECT_EQ(-45, ComputeEnergy(parsing::ParseDotBracketSecondary("GGGGAAACCCC", "((((...))))")).energy);
  EXPECT_EQ(72, ComputeEnergy(parsing::ParseDotBracketSecondary("UGACAAAGGCGA", "(..(...)...)")).energy);
}

}
}
