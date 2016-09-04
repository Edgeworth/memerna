#include <cstdlib>
#include "energy/load_model.h"
#include "fold/fold.h"
#include "parsing.h"
#include "gtest/gtest.h"
#include "common_test.h"

namespace memerna {
namespace fold {

class FoldAlgTest : public testing::TestWithParam<std::tuple<energy::EnergyModel, fold::context_options_t::TableAlg>> {
};


TEST_P(FoldAlgTest, T04) {
  const auto& em = std::get<0>(GetParam());
  fold::context_options_t options(std::get<1>(GetParam()));
  ONLY_FOR_THIS_MODEL(em, T04_MODEL_HASH);

  EXPECT_EQ(-45, fold::Context(parsing::StringToPrimary("GGGGAAACCCC"), em, options).Fold().energy);
  EXPECT_EQ(-51, fold::Context(parsing::StringToPrimary(
      "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"), em, options).Fold().energy);
  EXPECT_EQ(-133, fold::Context(parsing::StringToPrimary(
      "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"), em, options).Fold().energy);
  EXPECT_EQ(-57, fold::Context(parsing::StringToPrimary(
      "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"), em, options).Fold().energy);
  EXPECT_EQ(-121, fold::Context(parsing::StringToPrimary(
      "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"), em, options).Fold().energy);
  EXPECT_EQ(-74, fold::Context(parsing::StringToPrimary(
      "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"), em, options).Fold().energy);
  EXPECT_EQ(-32, fold::Context(parsing::StringToPrimary(
      "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"), em, options).Fold().energy);
  EXPECT_EQ(-40, fold::Context(parsing::StringToPrimary("CCGAAGGGGCUGCGGCG"), em, options).Fold().energy);
  EXPECT_EQ(-120, fold::Context(parsing::StringToPrimary(
      "CCGGGCCAGCCCGCUCCUACGGGGGGUC"), em, options).Fold().energy);
  EXPECT_EQ(-74, fold::Context(parsing::StringToPrimary(
      "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"), em, options).Fold().energy);
  EXPECT_EQ(-6, fold::Context(parsing::StringToPrimary("CCUCCGGG"), em, options).Fold().energy);
  EXPECT_EQ(-30, fold::Context(parsing::StringToPrimary(
      "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"), em, options).Fold().energy);
  EXPECT_EQ(-65, fold::Context(parsing::StringToPrimary(
      "CGCAGGGUCGGACCCGGGAGAACCGCGA"), em, options).Fold().energy);
  EXPECT_EQ(-60, fold::Context(parsing::StringToPrimary(
      "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"), em, options).Fold().energy);
  EXPECT_EQ(-6, fold::Context(parsing::StringToPrimary("CGGAAACGG"), em, options).Fold().energy);
  EXPECT_EQ(-22, fold::Context(parsing::StringToPrimary("CUGAAACUGGAAACAGAAAUG"), em, options).Fold().energy);
  EXPECT_EQ(-12, fold::Context(parsing::StringToPrimary("CUUAUAGUUAAGG"), em, options).Fold().energy);
  EXPECT_EQ(-122, fold::Context(parsing::StringToPrimary(
      "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"), em, options).Fold().energy);
  EXPECT_EQ(-29, fold::Context(parsing::StringToPrimary("GCCAAGGCCCCACCCGGA"), em, options).Fold().energy);
  EXPECT_EQ(-39, fold::Context(parsing::StringToPrimary(
      "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"), em, options).Fold().energy);
  EXPECT_EQ(-67, fold::Context(parsing::StringToPrimary(
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"), em, options).Fold().energy);
  EXPECT_EQ(-276, fold::Context(parsing::StringToPrimary(
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"), em, options).Fold().energy);
  EXPECT_EQ(-53, fold::Context(parsing::StringToPrimary(
      "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"), em, options).Fold().energy);
  EXPECT_EQ(-157, fold::Context(parsing::StringToPrimary(
      "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"), em, options).Fold().energy);
  EXPECT_EQ(-49, fold::Context(parsing::StringToPrimary("GGCCGAUGGCAGCGAUAGC"), em, options).Fold().energy);
  EXPECT_EQ(-44, fold::Context(parsing::StringToPrimary(
      "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"), em, options).Fold().energy);
  EXPECT_EQ(-45, fold::Context(parsing::StringToPrimary(
      "GGGGAAACCCC"), em, options).Fold().energy);
  EXPECT_EQ(-29, fold::Context(parsing::StringToPrimary(
      "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"), em, options).Fold().energy);
  EXPECT_EQ(-23, fold::Context(parsing::StringToPrimary(
      "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"), em, options).Fold().energy);
  EXPECT_EQ(-80, fold::Context(parsing::StringToPrimary(
      "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"), em, options).Fold().energy);
  EXPECT_EQ(-4, fold::Context(parsing::StringToPrimary(
      "UGCAAAGCAA"), em, options).Fold().energy);
  EXPECT_EQ(-208, fold::Context(parsing::StringToPrimary(
      "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"), em, options).Fold().energy);
}

INSTANTIATE_TEST_CASE_P(FoldAlgTest, FoldAlgTest,
    testing::Combine(
        testing::Values(energy::LoadEnergyModelFromDataDir(ENERGY_MODEL_PATH)),
        testing::ValuesIn(fold::context_options_t::TABLE_ALGS)
    ));

TEST(FoldTest, Precomp) {
  auto em = energy::LoadEnergyModelFromDataDir(ENERGY_MODEL_PATH);
  ONLY_FOR_THIS_MODEL(em, T04_MODEL_HASH);

  auto pc = internal::PrecomputeData(parsing::StringToPrimary("GGGGAAACCCC"), em);
  EXPECT_EQ(-21 - 4 - 16, pc.min_mismatch_coax);
  EXPECT_EQ(-34, pc.min_flush_coax);
  EXPECT_EQ(-26, pc.min_twoloop_not_stack);

  energy_t augubranch[4][4] = {
      {-6, -6, -6, 5 - 6},
      {-6, -6, -6, -6},
      {-6, -6, -6, 5 - 6},
      {5 - 6, -6, 5 - 6, -6}
  };
  EXPECT_EQ(sizeof(augubranch), sizeof(pc.augubranch));
  EXPECT_TRUE(std::memcmp(augubranch, pc.augubranch, sizeof(augubranch)) == 0);
}

TEST(FoldTest, Helpers) {
  EXPECT_EQ(0, internal::MaxNumContiguous(parsing::StringToPrimary("")));
  EXPECT_EQ(1, internal::MaxNumContiguous(parsing::StringToPrimary("A")));
  EXPECT_EQ(2, internal::MaxNumContiguous(parsing::StringToPrimary("AA")));
  EXPECT_EQ(2, internal::MaxNumContiguous(parsing::StringToPrimary("GUAAC")));
  EXPECT_EQ(1, internal::MaxNumContiguous(parsing::StringToPrimary("GUACA")));
  EXPECT_EQ(3, internal::MaxNumContiguous(parsing::StringToPrimary("GAUCCC")));
  EXPECT_EQ(3, internal::MaxNumContiguous(parsing::StringToPrimary("GGGAUC")));
  EXPECT_EQ(4, internal::MaxNumContiguous(parsing::StringToPrimary("GGGAUCAAAA")));
  EXPECT_EQ(5, internal::MaxNumContiguous(parsing::StringToPrimary("GGGAUUUUUCAAAA")));
}

}
}
