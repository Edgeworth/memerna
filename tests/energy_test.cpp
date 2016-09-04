#include "constants.h"
#include "fold/fold.h"
#include "parsing.h"
#include "gtest/gtest.h"
#include "common_test.h"
#include "energy/load_model.h"

namespace memerna {
namespace energy {

class EnergyTest : public testing::TestWithParam<energy::EnergyModel> {
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


TEST_P(EnergyTest, MultiloopEnergy) {
  const auto& em = GetParam();
  EXPECT_EQ(em.multiloop_hack_a + 4 * em.multiloop_hack_b, em.MultiloopInitiation(4));
}

TEST_P(EnergyTest, NNDBHairpinLoopExamples) {
  const auto& em = GetParam();
  EXPECT_EQ(em.stack[C][A][U][G] + em.stack[A][C][G][U] + em.stack[C][A][U][G] + em.augu_penalty +
      em.terminal[A][A][A][U] + em.HairpinInitiation(6),
      ComputeEnergy(kNNDBHairpin1, em).energy);
  EXPECT_EQ(em.stack[C][A][U][G] + em.stack[A][C][G][U] + em.stack[C][A][U][G] + em.augu_penalty +
      em.terminal[A][G][G][U] + em.hairpin_gg_first_mismatch + em.HairpinInitiation(5),
      ComputeEnergy(kNNDBHairpin2, em).energy);

  auto iter = em.hairpin.find("CCGAGG");
  auto hairpin_val = (iter == em.hairpin.end() ? 0 : iter->second);
  EXPECT_EQ(em.stack[C][A][U][G] + em.stack[A][C][G][U] + em.stack[C][C][G][G] + hairpin_val,
      ComputeEnergy(kNNDBHairpin3, em).energy);
  EXPECT_EQ(em.stack[C][A][U][G] + em.stack[A][C][G][U] + em.stack[C][A][U][G] +
      em.augu_penalty + em.terminal[A][C][C][U] + em.HairpinInitiation(6) +
      em.hairpin_all_c_a * 6 + em.hairpin_all_c_b, ComputeEnergy(kNNDBHairpin4, em).energy);
  EXPECT_EQ(em.stack[C][G][C][G] + em.stack[G][G][C][C] + em.stack[G][G][U][C] + em.augu_penalty +
      em.terminal[G][G][G][U] + em.hairpin_gg_first_mismatch + em.HairpinInitiation(5) + em.hairpin_special_gu_closure,
      ComputeEnergy(kNNDBHairpin5, em).energy);

  EXPECT_EQ(em.augu_penalty + em.terminal[A][A][A][U] + em.HairpinInitiation(6),
      fold::Context(kNNDBHairpin1.r, em).FastHairpin(3, 10));
  EXPECT_EQ(
      em.augu_penalty + em.terminal[A][G][G][U] + em.hairpin_gg_first_mismatch + em.HairpinInitiation(5),
      fold::Context(kNNDBHairpin2.r, em).FastHairpin(3, 9));
  EXPECT_EQ(hairpin_val,
      fold::Context(kNNDBHairpin3.r, em).FastHairpin(3, 8));
  EXPECT_EQ(em.augu_penalty + em.terminal[A][C][C][U] + em.HairpinInitiation(6) +
      em.hairpin_all_c_a * 6 + em.hairpin_all_c_b,
      fold::Context(kNNDBHairpin4.r, em).FastHairpin(3, 10));
  EXPECT_EQ(em.augu_penalty + em.terminal[G][G][G][U] + em.hairpin_gg_first_mismatch +
      em.HairpinInitiation(5) + em.hairpin_special_gu_closure,
      fold::Context(kNNDBHairpin5.r, em).FastHairpin(3, 9));
}

TEST_P(EnergyTest, NNDBBulgeLoopExamples) {
  const auto& em = GetParam();
  EXPECT_EQ(em.stack[G][C][G][C] + em.stack[C][C][G][G] + em.BulgeInitiation(1) +
      em.bulge_special_c + em.stack[C][G][C][G] + em.HairpinInitiation(3) -
      energy_t(round(10.0 * constants::R * constants::T * log(3))),
      ComputeEnergy(kNNDBBulge1, em).energy);

  EXPECT_EQ(em.stack[G][A][U][C] + em.augu_penalty + em.BulgeInitiation(3) + em.HairpinInitiation(3),
      ComputeEnergy(kNNDBBulge2, em).energy);
}


TEST_P(EnergyTest, NNDBMultiloopExamples) {
  const auto& em = GetParam();
  EXPECT_EQ(em.stack[C][A][U][G] + em.stack[A][C][G][U] + em.stack[C][A][U][G] +
      2 * em.augu_penalty + 2 * em.HairpinInitiation(3),
      ComputeEnergy(kFlushCoax, em).energy);
  EXPECT_EQ(em.stack[G][A][U][C] + em.terminal[C][G][A][G] + em.coax_mismatch_non_contiguous +
      3 * em.HairpinInitiation(3) + em.MultiloopInitiation(4) + 2 * em.augu_penalty,
      ComputeEnergy(kNNDBMultiloop, em).energy);
}

TEST_P(EnergyTest, NNDBInternalLoopExamples) {
  const auto& em = GetParam();
  EXPECT_EQ(em.stack[C][A][U][G] + em.stack[C][G][C][G] + em.InternalLoopInitiation(5) +
      std::min(em.internal_asym, constants::NINIO_MAX_ASYM) + em.internal_2x3_mismatch[A][G][G][U] +
      em.internal_2x3_mismatch[G][G][A][C] + em.internal_augu_penalty + em.HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x3, em).energy);
  EXPECT_EQ(em.stack[C][A][U][G] + em.stack[C][G][C][G] + em.internal_2x2[A][G][A][C][G][G][A][U] +
      em.HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal2x2, em).energy);
  EXPECT_EQ(em.stack[C][A][U][G] + em.stack[C][G][C][G] + em.InternalLoopInitiation(6) +
      std::min(4 * em.internal_asym, constants::NINIO_MAX_ASYM) +
      em.internal_augu_penalty + em.HairpinInitiation(3),
      ComputeEnergy(kNNDBInternal1x5, em).energy);
}

TEST_P(EnergyTest, BaseCases) {
  const auto& em = GetParam();
  EXPECT_EQ(em.augu_penalty + em.stack[G][A][U][C] + em.hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketSecondary("GAAAAUC", "((...))"), em).energy);
  EXPECT_EQ(em.augu_penalty * 2 + em.stack[G][A][U][U] + em.hairpin_init[3],
      ComputeEnergy(parsing::ParseDotBracketSecondary("GAAAAUU", "((...))"), em).energy);
  EXPECT_EQ(em.augu_penalty * 2 + em.HairpinInitiation(3) +
      std::min(em.terminal[U][A][A][A], std::min(em.dangle3[U][A][A], em.dangle5[U][A][A])),
      ComputeEnergy(parsing::ParseDotBracketSecondary("AAAAAUA", ".(...)."), em).energy);
  EXPECT_EQ(em.augu_penalty * 2 + em.HairpinInitiation(3),
      ComputeEnergy(parsing::ParseDotBracketSecondary("AAAAU", "(...)"), em).energy);
  EXPECT_EQ(em.stack[G][C][G][C] + em.stack[C][U][A][G] + em.BulgeInitiation(1) +
      em.stack[U][G][C][A] + em.HairpinInitiation(3),
      ComputeEnergy(kBulge1, em).energy);
  EXPECT_EQ(em.InternalLoopInitiation(5) + em.internal_asym + em.internal_augu_penalty + em.augu_penalty +
      em.internal_2x3_mismatch[A][G][A][U] + em.internal_2x3_mismatch[C][A][A][G] + em.HairpinInitiation(3),
      ComputeEnergy(kInternal1, em).energy);
}

TEST_P(EnergyTest, T04Tests) {
  const auto& em = GetParam();
  ONLY_FOR_THIS_MODEL(em, T04_MODEL_HASH);

  EXPECT_EQ(88, em.HairpinInitiation(87));
  EXPECT_EQ(68, em.BulgeInitiation(57));
  EXPECT_EQ(46, em.InternalLoopInitiation(67));

  EXPECT_EQ(45, ComputeEnergy(parsing::ParseDotBracketSecondary("GCAAAGCC", "((...).)"), em).energy);
  EXPECT_EQ(57, ComputeEnergy(parsing::ParseDotBracketSecondary("CCCAAAAUG", ".(.(...))"), em).energy);
  EXPECT_EQ(55, ComputeEnergy(parsing::ParseDotBracketSecondary("UACAGA", "(....)"), em).energy);
  EXPECT_EQ(-6, ComputeEnergy(parsing::ParseDotBracketSecondary("AGGGUCAUCCG", ".(((...)))."), em).energy);
  EXPECT_EQ(80, ComputeEnergy(parsing::ParseDotBracketSecondary("AGAGAAACAAAU", "(..(...)...)"), em).energy);
  EXPECT_EQ(95, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"), em).energy);
  EXPECT_EQ(77, ComputeEnergy(parsing::ParseDotBracketSecondary("CCCGAAACAG", "(..(...).)"), em).energy);
  EXPECT_EQ(74, ComputeEnergy(parsing::ParseDotBracketSecondary("GACAGAAACGCUGAAUC", "((..(...)......))"), em).energy);
  EXPECT_EQ(173, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"), em).energy);
  EXPECT_EQ(182, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"), em).energy);
  EXPECT_EQ(176, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"), em).energy);
  EXPECT_EQ(131, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"), em).energy);
  EXPECT_EQ(-276, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
      "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."), em).energy);
  EXPECT_EQ(179, ComputeEnergy(parsing::ParseDotBracketSecondary(
      "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"), em).energy);

  // Special stacking - this is not implemented. TODO: Implement this?
  EXPECT_EQ(37, ComputeEnergy(parsing::ParseDotBracketSecondary("GGUCAAAGGUC", "((((...))))"), em).energy);
  EXPECT_EQ(-45, ComputeEnergy(parsing::ParseDotBracketSecondary("GGGGAAACCCC", "((((...))))"), em).energy);
  EXPECT_EQ(72, ComputeEnergy(parsing::ParseDotBracketSecondary("UGACAAAGGCGA", "(..(...)...)"), em).energy);
}

INSTANTIATE_TEST_CASE_P(EnergyTest, EnergyTest, testing::Values(
    energy::LoadEnergyModelFromDataDir(ENERGY_MODEL_PATH)
));

}
}
