#include "constants.h"
#include "gtest/gtest.h"
#include "common_test.h"
#include "fold/context.h"
#include "fold/globals.h"
#include "parsing.h"

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
  EXPECT_EQ(g_em->multiloop_hack_a + 4 * g_em->multiloop_hack_b, g_em->MultiloopInitiation(4));
}

TEST_F(EnergyTest, NNDBHairpinLoopExamples) {
  EXPECT_EQ(g_em->stack[C][A][U][G] + g_em->stack[A][C][G][U] + g_em->stack[C][A][U][G] + g_em->augu_penalty +
      g_em->terminal[A][A][A][U] + g_em->HairpinInitiation(6),
      ComputeEnergy(kNNDBHairpin1, *g_em).energy);
  EXPECT_EQ(g_em->stack[C][A][U][G] + g_em->stack[A][C][G][U] + g_em->stack[C][A][U][G] + g_em->augu_penalty +
      g_em->terminal[A][G][G][U] + g_em->hairpin_gg_first_mismatch + g_em->HairpinInitiation(5),
      ComputeEnergy(kNNDBHairpin2, *g_em).energy);

  EXPECT_EQ(g_em->stack[C][A][U][G] + g_em->stack[A][C][G][U] + g_em->stack[C][C][G][G] + g_em->hairpin["CCGAGG"],
      ComputeEnergy(kNNDBHairpin3, *g_em).energy);
  EXPECT_EQ(g_em->stack[C][A][U][G] + g_em->stack[A][C][G][U] + g_em->stack[C][A][U][G] +
      g_em->augu_penalty + g_em->terminal[A][C][C][U] + g_em->HairpinInitiation(6) +
      g_em->hairpin_all_c_a * 6 + g_em->hairpin_all_c_b, ComputeEnergy(kNNDBHairpin4, *g_em).energy);
  EXPECT_EQ(g_em->stack[C][G][C][G] + g_em->stack[G][G][C][C] + g_em->stack[G][G][U][C] + g_em->augu_penalty +
      g_em->terminal[G][G][G][U] + g_em->hairpin_gg_first_mismatch + g_em->HairpinInitiation(5) +
      g_em->hairpin_special_gu_closure,
      ComputeEnergy(kNNDBHairpin5, *g_em).energy);

  fold::internal::SetGlobalState(kNNDBHairpin1.r, *g_em);
  EXPECT_EQ(g_em->augu_penalty + g_em->terminal[A][A][A][U] + g_em->HairpinInitiation(6),
      fold::internal::FastHairpin(3, 10));
  fold::internal::SetGlobalState(kNNDBHairpin2.r, *g_em);
  EXPECT_EQ(
      g_em->augu_penalty + g_em->terminal[A][G][G][U] + g_em->hairpin_gg_first_mismatch + g_em->HairpinInitiation(5),
      fold::internal::FastHairpin(3, 9));
  fold::internal::SetGlobalState(kNNDBHairpin3.r, *g_em);
  EXPECT_EQ(g_em->hairpin["CCGAGG"],
      fold::internal::FastHairpin(3, 8));
  fold::internal::SetGlobalState(kNNDBHairpin4.r, *g_em);
  EXPECT_EQ(g_em->augu_penalty + g_em->terminal[A][C][C][U] + g_em->HairpinInitiation(6) +
      g_em->hairpin_all_c_a * 6 + g_em->hairpin_all_c_b,
      fold::internal::FastHairpin(3, 10));
  fold::internal::SetGlobalState(kNNDBHairpin5.r, *g_em);
  EXPECT_EQ(g_em->augu_penalty + g_em->terminal[G][G][G][U] + g_em->hairpin_gg_first_mismatch +
      g_em->HairpinInitiation(5) + g_em->hairpin_special_gu_closure,
      fold::internal::FastHairpin(3, 9));
}

TEST_F(EnergyTest, NNDBBulgeLoopExamples) {
  EXPECT_EQ(g_em->stack[G][C][G][C] + g_em->stack[C][C][G][G] + g_em->BulgeInitiation(1) +
      g_em->bulge_special_c + g_em->stack[C][G][C][G] + g_em->HairpinInitiation(3) -
      energy_t(round(10.0 * constants::R * constants::T * log(3))),
      ComputeEnergy(kNNDBBulge1, *g_em).energy);

  EXPECT_EQ(g_em->stack[G][A][U][C] + g_em->augu_penalty + g_em->BulgeInitiation(3) + g_em->HairpinInitiation(3),
      ComputeEnergy(kNNDBBulge2, *g_em).energy);
}


TEST_F(EnergyTest, NNDBMultiloopExamples) {
  EXPECT_EQ(g_em->stack[C][A][U][G] + g_em->stack[A][C][G][U] + g_em->stack[C][A][U][G] +
      2 * g_em->augu_penalty + 2 * g_em->HairpinInitiation(3),
      ComputeEnergy(kFlushCoax, *g_em).energy);
  EXPECT_EQ(g_em->stack[G][A][U][C] + g_em->terminal[C][G][A][G] + g_em->coax_mismatch_non_contiguous +
      3 * g_em->HairpinInitiation(3) + g_em->MultiloopInitiation(4) + 2 * g_em->augu_penalty,
      ComputeEnergy(kNNDBMultiloop, *g_em).energy);
}

TEST_F(EnergyTest, NNDBInternalLoopExamples) {
  EXPECT_EQ(g_em->stack[C][A][U][G] + g_em->stack[C][G][C][G] + g_em->InternalLoopInitiation(5) +
      std::min(g_em->internal_asym, constants::NINIO_MAX_ASYM) + g_em->internal_2x3_mismatch[A][G][G][U] +
      g_em->internal_2x3_mismatch[G][G][A][C] + g_em->internal_augu_penalty + g_em->HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x3, *g_em).energy);
  EXPECT_EQ(g_em->stack[C][A][U][G] + g_em->stack[C][G][C][G] + g_em->internal_2x2[A][G][A][C][G][G][A][U] +
      g_em->HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x2, *g_em).energy);
  EXPECT_EQ(g_em->stack[C][A][U][G] + g_em->stack[C][G][C][G] + g_em->InternalLoopInitiation(6) +
      std::min(4 * g_em->internal_asym, constants::NINIO_MAX_ASYM) +
      g_em->internal_augu_penalty + g_em->HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal1x5, *g_em).energy);
}

TEST_F(EnergyTest, BaseCases) {
  EXPECT_EQ(g_em->augu_penalty + g_em->stack[G][A][U][C] + g_em->hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketSecondary("GAAAAUC", "((...))"), *g_em).energy);
  EXPECT_EQ(g_em->augu_penalty * 2 + g_em->stack[G][A][U][U] + g_em->hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketSecondary("GAAAAUU", "((...))"), *g_em).energy);
  EXPECT_EQ(g_em->augu_penalty * 2 + g_em->HairpinInitiation(3) +
      std::min(g_em->terminal[U][A][A][A], std::min(g_em->dangle3[U][A][A], g_em->dangle5[U][A][A])),
      ComputeEnergy(parsing::ParseDotBracketSecondary("AAAAAUA", ".(...)."), *g_em).energy);
  EXPECT_EQ(g_em->augu_penalty * 2 + g_em->HairpinInitiation(3),
      ComputeEnergy(parsing::ParseDotBracketSecondary("AAAAU", "(...)"), *g_em).energy);
  EXPECT_EQ(g_em->stack[G][C][G][C] + g_em->stack[C][U][A][G] + g_em->BulgeInitiation(1) +
      g_em->stack[U][G][C][A] + g_em->HairpinInitiation(3),
      ComputeEnergy(kBulge1, *g_em).energy);
  EXPECT_EQ(g_em->InternalLoopInitiation(5) + g_em->internal_asym + g_em->internal_augu_penalty + g_em->augu_penalty +
      g_em->internal_2x3_mismatch[A][G][A][U] + g_em->internal_2x3_mismatch[C][A][A][G] + g_em->HairpinInitiation(3),
      ComputeEnergy(kInternal1, *g_em).energy);
}

TEST_F(EnergyTest, T04Tests) {
  ONLY_FOR_THIS_MODEL(g_em, T04_MODEL_HASH);

  EXPECT_EQ(88, g_em->HairpinInitiation(87));
  EXPECT_EQ(68, g_em->BulgeInitiation(57));
  EXPECT_EQ(46, g_em->InternalLoopInitiation(67));

  EXPECT_EQ(45, ComputeEnergy(parsing::ParseDotBracketSecondary("GCAAAGCC", "((...).)"), *g_em).energy);
  EXPECT_EQ(57, ComputeEnergy(parsing::ParseDotBracketSecondary("CCCAAAAUG", ".(.(...))"), *g_em).energy);
  EXPECT_EQ(55, ComputeEnergy(parsing::ParseDotBracketSecondary("UACAGA", "(....)"), *g_em).energy);
  EXPECT_EQ(-6, ComputeEnergy(parsing::ParseDotBracketSecondary("AGGGUCAUCCG", ".(((...)))."), *g_em).energy);
  EXPECT_EQ(80, ComputeEnergy(parsing::ParseDotBracketSecondary("AGAGAAACAAAU", "(..(...)...)"), *g_em).energy);
  EXPECT_EQ(95, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"), *g_em).energy);
  EXPECT_EQ(77, ComputeEnergy(parsing::ParseDotBracketSecondary("CCCGAAACAG", "(..(...).)"), *g_em).energy);
  EXPECT_EQ(74,
      ComputeEnergy(parsing::ParseDotBracketSecondary("GACAGAAACGCUGAAUC", "((..(...)......))"), *g_em).energy);
  EXPECT_EQ(173, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"), *g_em).energy);
  EXPECT_EQ(182, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"), *g_em).energy);
  EXPECT_EQ(176, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"), *g_em).energy);
  EXPECT_EQ(131, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"), *g_em).energy);
  EXPECT_EQ(-276, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
      "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."), *g_em).energy);
  EXPECT_EQ(179, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"), *g_em).energy);

  // Special stacking - this is not implg_emented. TODO: Implg_ement this?
  EXPECT_EQ(37, ComputeEnergy(parsing::ParseDotBracketSecondary("GGUCAAAGGUC", "((((...))))"), *g_em).energy);
  EXPECT_EQ(-45, ComputeEnergy(parsing::ParseDotBracketSecondary("GGGGAAACCCC", "((((...))))"), *g_em).energy);
  EXPECT_EQ(72, ComputeEnergy(parsing::ParseDotBracketSecondary("UGACAAAGGCGA", "(..(...)...)"), *g_em).energy);
}

}
}
