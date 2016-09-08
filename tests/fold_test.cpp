#include <cstdlib>
#include "common_test.h"
#include "gtest/gtest.h"
#include "fold/context.h"
#include "fold/precomp.h"
#include "parsing.h"

namespace memerna {
namespace fold {

class FoldAlgTest : public testing::TestWithParam<fold::context_options_t::TableAlg> {
};


TEST_P(FoldAlgTest, T04) {
  fold::context_options_t options(GetParam());
  ONLY_FOR_THIS_MODEL(g_em, T04_MODEL_HASH);

  EXPECT_EQ(-45, fold::Context(parsing::StringToPrimary("GGGGAAACCCC"), g_em, options).Fold().energy);
  EXPECT_EQ(-51, fold::Context(parsing::StringToPrimary(
      "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"), g_em, options).Fold().energy);
  EXPECT_EQ(-133, fold::Context(parsing::StringToPrimary(
      "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"), g_em, options).Fold().energy);
  EXPECT_EQ(-57, fold::Context(parsing::StringToPrimary(
      "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"), g_em, options).Fold().energy);
  EXPECT_EQ(-121, fold::Context(parsing::StringToPrimary(
      "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"), g_em, options).Fold().energy);
  EXPECT_EQ(-74, fold::Context(parsing::StringToPrimary(
      "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"), g_em, options).Fold().energy);
  EXPECT_EQ(-32, fold::Context(parsing::StringToPrimary(
      "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"), g_em, options).Fold().energy);
  EXPECT_EQ(-40, fold::Context(parsing::StringToPrimary("CCGAAGGGGCUGCGGCG"), g_em, options).Fold().energy);
  EXPECT_EQ(-120, fold::Context(parsing::StringToPrimary(
      "CCGGGCCAGCCCGCUCCUACGGGGGGUC"), g_em, options).Fold().energy);
  EXPECT_EQ(-74, fold::Context(parsing::StringToPrimary(
      "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"), g_em, options).Fold().energy);
  EXPECT_EQ(-6, fold::Context(parsing::StringToPrimary("CCUCCGGG"), g_em, options).Fold().energy);
  EXPECT_EQ(-30, fold::Context(parsing::StringToPrimary(
      "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"), g_em, options).Fold().energy);
  EXPECT_EQ(-65, fold::Context(parsing::StringToPrimary(
      "CGCAGGGUCGGACCCGGGAGAACCGCGA"), g_em, options).Fold().energy);
  EXPECT_EQ(-60, fold::Context(parsing::StringToPrimary(
      "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"), g_em, options).Fold().energy);
  EXPECT_EQ(-6, fold::Context(parsing::StringToPrimary("CGGAAACGG"), g_em, options).Fold().energy);
  EXPECT_EQ(-22, fold::Context(parsing::StringToPrimary("CUGAAACUGGAAACAGAAAUG"), g_em, options).Fold().energy);
  EXPECT_EQ(-12, fold::Context(parsing::StringToPrimary("CUUAUAGUUAAGG"), g_em, options).Fold().energy);
  EXPECT_EQ(-122, fold::Context(parsing::StringToPrimary(
      "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"), g_em, options).Fold().energy);
  EXPECT_EQ(-29, fold::Context(parsing::StringToPrimary("GCCAAGGCCCCACCCGGA"), g_em, options).Fold().energy);
  EXPECT_EQ(-39, fold::Context(parsing::StringToPrimary(
      "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"), g_em, options).Fold().energy);
  EXPECT_EQ(-67, fold::Context(parsing::StringToPrimary(
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"), g_em, options).Fold().energy);
  EXPECT_EQ(-276, fold::Context(parsing::StringToPrimary(
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"), g_em, options).Fold().energy);
  EXPECT_EQ(-53, fold::Context(parsing::StringToPrimary(
      "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"), g_em, options).Fold().energy);
  EXPECT_EQ(-157, fold::Context(parsing::StringToPrimary(
      "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"), g_em, options).Fold().energy);
  EXPECT_EQ(-49, fold::Context(parsing::StringToPrimary("GGCCGAUGGCAGCGAUAGC"), g_em, options).Fold().energy);
  EXPECT_EQ(-44, fold::Context(parsing::StringToPrimary(
      "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"), g_em, options).Fold().energy);
  EXPECT_EQ(-45, fold::Context(parsing::StringToPrimary(
      "GGGGAAACCCC"), g_em, options).Fold().energy);
  EXPECT_EQ(-29, fold::Context(parsing::StringToPrimary(
      "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"), g_em, options).Fold().energy);
  EXPECT_EQ(-23, fold::Context(parsing::StringToPrimary(
      "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"), g_em, options).Fold().energy);
  EXPECT_EQ(-80, fold::Context(parsing::StringToPrimary(
      "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"), g_em, options).Fold().energy);
  EXPECT_EQ(-4, fold::Context(parsing::StringToPrimary(
      "UGCAAAGCAA"), g_em, options).Fold().energy);
  EXPECT_EQ(-208, fold::Context(parsing::StringToPrimary(
      "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"), g_em, options).Fold().energy);
}

INSTANTIATE_TEST_CASE_P(FoldAlgTest, FoldAlgTest, testing::ValuesIn(fold::context_options_t::TABLE_ALGS));

TEST(FoldTest, Precomp) {
  ONLY_FOR_THIS_MODEL(g_em, T04_MODEL_HASH);

  auto pc = internal::PrecomputeData(parsing::StringToPrimary("GGGGAAACCCC"), *g_em);
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
