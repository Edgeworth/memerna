#include <cstdlib>
#include "fold/fold.h"
#include "fold/fold_globals.h"
#include "parsing.h"
#include "gtest/gtest.h"
#include "common_test.h"
#include "energy/energy_model.h"

namespace memerna {
namespace fold {

class FoldAlgTest : public testing::TestWithParam<fold_fn_t*> {
};


TEST_P(FoldAlgTest, T04) {
  ONLY_FOR_THIS_MODEL(T04_MODEL_HASH);

  auto fold_fn = [this](const auto& primary) {return GetParam()(primary, nullptr);};
  EXPECT_EQ(-45, fold_fn(parsing::StringToPrimary("GGGGAAACCCC")).energy);
  EXPECT_EQ(-51, fold_fn(parsing::StringToPrimary(
      "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA")).energy);
  EXPECT_EQ(-133, fold_fn(parsing::StringToPrimary(
      "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU")).energy);
  EXPECT_EQ(-57, fold_fn(parsing::StringToPrimary("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC")).energy);
  EXPECT_EQ(-121, fold_fn(parsing::StringToPrimary("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG")).energy);
  EXPECT_EQ(-74, fold_fn(parsing::StringToPrimary("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG")).energy);
  EXPECT_EQ(-32, fold_fn(parsing::StringToPrimary("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG")).energy);
  EXPECT_EQ(-40, fold_fn(parsing::StringToPrimary("CCGAAGGGGCUGCGGCG")).energy);
  EXPECT_EQ(-120, fold_fn(parsing::StringToPrimary("CCGGGCCAGCCCGCUCCUACGGGGGGUC")).energy);
  EXPECT_EQ(-74, fold_fn(parsing::StringToPrimary("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG")).energy);
  EXPECT_EQ(-6, fold_fn(parsing::StringToPrimary("CCUCCGGG")).energy);
  EXPECT_EQ(-30, fold_fn(parsing::StringToPrimary("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC")).energy);
  EXPECT_EQ(-65, fold_fn(parsing::StringToPrimary("CGCAGGGUCGGACCCGGGAGAACCGCGA")).energy);
  EXPECT_EQ(-60, fold_fn(parsing::StringToPrimary("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA")).energy);
  EXPECT_EQ(-6, fold_fn(parsing::StringToPrimary("CGGAAACGG")).energy);
  EXPECT_EQ(-22, fold_fn(parsing::StringToPrimary("CUGAAACUGGAAACAGAAAUG")).energy);
  EXPECT_EQ(-12, fold_fn(parsing::StringToPrimary("CUUAUAGUUAAGG")).energy);
  EXPECT_EQ(-122, fold_fn(parsing::StringToPrimary(
      "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC")).energy);
  EXPECT_EQ(-29, fold_fn(parsing::StringToPrimary("GCCAAGGCCCCACCCGGA")).energy);
  EXPECT_EQ(-39, fold_fn(parsing::StringToPrimary("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC")).energy);
  EXPECT_EQ(-67, fold_fn(parsing::StringToPrimary("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC")).energy);
  EXPECT_EQ(-276, fold_fn(parsing::StringToPrimary(
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA")).energy);
  EXPECT_EQ(-53, fold_fn(parsing::StringToPrimary("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU")).energy);
  EXPECT_EQ(-157, fold_fn(parsing::StringToPrimary("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG")).energy);
  EXPECT_EQ(-49, fold_fn(parsing::StringToPrimary("GGCCGAUGGCAGCGAUAGC")).energy);
  EXPECT_EQ(-44, fold_fn(parsing::StringToPrimary("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC")).energy);
  EXPECT_EQ(-45, fold_fn(parsing::StringToPrimary("GGGGAAACCCC")).energy);
  EXPECT_EQ(-29, fold_fn(parsing::StringToPrimary("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG")).energy);
  EXPECT_EQ(-23, fold_fn(parsing::StringToPrimary("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG")).energy);
  EXPECT_EQ(-80, fold_fn(parsing::StringToPrimary("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG")).energy);
  EXPECT_EQ(-4, fold_fn(parsing::StringToPrimary("UGCAAAGCAA")).energy);
  EXPECT_EQ(-208, fold_fn(parsing::StringToPrimary(
      "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU")).energy);
}

INSTANTIATE_TEST_CASE_P(FoldAlgTest, FoldAlgTest, testing::ValuesIn(FOLD_FUNCTIONS));

TEST(FoldTest, Constants) {
  ONLY_FOR_THIS_MODEL(T04_MODEL_HASH);

  r = parsing::StringToPrimary("GGGGAAACCCC");
  fold::InitFold();
  EXPECT_EQ(-21 - 4 - 16, g_min_mismatch_coax);
  EXPECT_EQ(-34, g_min_flush_coax);
  EXPECT_EQ(-26, g_min_twoloop_not_stack);

  energy_t augubranch[4][4] = {
      {-6, -6, -6, 5 - 6},
      {-6, -6, -6, -6},
      {-6, -6, -6, 5 - 6},
      {5 - 6, -6, 5 - 6, -6}
  };
  EXPECT_EQ(sizeof(augubranch), sizeof(g_augubranch));
  EXPECT_TRUE(std::memcmp(augubranch, g_augubranch, sizeof(augubranch)) == 0);
}

TEST(FoldTest, Helpers) {
  EXPECT_EQ(0, fold::MaxNumContiguous(parsing::StringToPrimary("")));
  EXPECT_EQ(1, fold::MaxNumContiguous(parsing::StringToPrimary("A")));
  EXPECT_EQ(2, fold::MaxNumContiguous(parsing::StringToPrimary("AA")));
  EXPECT_EQ(2, fold::MaxNumContiguous(parsing::StringToPrimary("GUAAC")));
  EXPECT_EQ(1, fold::MaxNumContiguous(parsing::StringToPrimary("GUACA")));
  EXPECT_EQ(3, fold::MaxNumContiguous(parsing::StringToPrimary("GAUCCC")));
  EXPECT_EQ(3, fold::MaxNumContiguous(parsing::StringToPrimary("GGGAUC")));
  EXPECT_EQ(4, fold::MaxNumContiguous(parsing::StringToPrimary("GGGAUCAAAA")));
  EXPECT_EQ(5, fold::MaxNumContiguous(parsing::StringToPrimary("GGGAUUUUUCAAAA")));
}

}
}
