#include <cstdlib>
#include "fold/fold.h"
#include "fold/fold_globals.h"
#include "parsing.h"
#include "gtest/gtest.h"
#include "common_test.h"

namespace memerna {
namespace fold {

class FoldAlgTest : public testing::TestWithParam<fold_fn_t*> {
};


TEST_P(FoldAlgTest, T04) {
  if (EnergyModelChecksum() != T04_MODEL_HASH) {
    printf("Skipping energy model specific energy tests.");
    return;
  }

  auto fold_fn = [this](const auto& rna) {return GetParam()(rna, nullptr);};

  EXPECT_EQ(-45, fold_fn(parsing::StringToRna("GGGGAAACCCC")).energy);
  EXPECT_EQ(-51, fold_fn(parsing::StringToRna("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA")).energy);
  EXPECT_EQ(-133,
      fold_fn(parsing::StringToRna("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU")).energy);
  EXPECT_EQ(-57, fold_fn(parsing::StringToRna("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC")).energy);
  EXPECT_EQ(-121, fold_fn(parsing::StringToRna("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG")).energy);
  EXPECT_EQ(-74, fold_fn(parsing::StringToRna("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG")).energy);
  EXPECT_EQ(-32, fold_fn(parsing::StringToRna("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG")).energy);
  EXPECT_EQ(-40, fold_fn(parsing::StringToRna("CCGAAGGGGCUGCGGCG")).energy);
  EXPECT_EQ(-120, fold_fn(parsing::StringToRna("CCGGGCCAGCCCGCUCCUACGGGGGGUC")).energy);
  EXPECT_EQ(-74, fold_fn(parsing::StringToRna("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG")).energy);
  EXPECT_EQ(-6, fold_fn(parsing::StringToRna("CCUCCGGG")).energy);
  EXPECT_EQ(-30, fold_fn(parsing::StringToRna("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC")).energy);
  EXPECT_EQ(-65, fold_fn(parsing::StringToRna("CGCAGGGUCGGACCCGGGAGAACCGCGA")).energy);
  EXPECT_EQ(-60, fold_fn(parsing::StringToRna("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA")).energy);
  EXPECT_EQ(-6, fold_fn(parsing::StringToRna("CGGAAACGG")).energy);
  EXPECT_EQ(-22, fold_fn(parsing::StringToRna("CUGAAACUGGAAACAGAAAUG")).energy);
  EXPECT_EQ(-12, fold_fn(parsing::StringToRna("CUUAUAGUUAAGG")).energy);
  EXPECT_EQ(-122,
      fold_fn(parsing::StringToRna("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC")).energy);
  EXPECT_EQ(-29, fold_fn(parsing::StringToRna("GCCAAGGCCCCACCCGGA")).energy);
  EXPECT_EQ(-39, fold_fn(parsing::StringToRna("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC")).energy);
  EXPECT_EQ(-67, fold_fn(parsing::StringToRna("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC")).energy);
  EXPECT_EQ(-276,
      fold_fn(parsing::StringToRna(
          "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA")).energy);
  EXPECT_EQ(-53, fold_fn(parsing::StringToRna("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU")).energy);
  EXPECT_EQ(-157, fold_fn(parsing::StringToRna("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG")).energy);
  EXPECT_EQ(-49, fold_fn(parsing::StringToRna("GGCCGAUGGCAGCGAUAGC")).energy);
  EXPECT_EQ(-44, fold_fn(parsing::StringToRna("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC")).energy);
  EXPECT_EQ(-45, fold_fn(parsing::StringToRna("GGGGAAACCCC")).energy);
  EXPECT_EQ(-29, fold_fn(parsing::StringToRna("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG")).energy);
  EXPECT_EQ(-23, fold_fn(parsing::StringToRna("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG")).energy);
  EXPECT_EQ(-80, fold_fn(parsing::StringToRna("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG")).energy);
  EXPECT_EQ(-4, fold_fn(parsing::StringToRna("UGCAAAGCAA")).energy);
  EXPECT_EQ(-208, fold_fn(parsing::StringToRna("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU")).energy);
}

INSTANTIATE_TEST_CASE_P(FoldAlgTest, FoldAlgTest, testing::Values(&Fold0, &Fold1, &Fold2, &Fold3));

TEST(FoldTest, Constants) {
  if (EnergyModelChecksum() != T04_MODEL_HASH) {
    printf("Skipping energy model specific energy tests.");
    return;
  }

  r = parsing::StringToRna("GGGGAAACCCC");
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
  EXPECT_EQ(0, fold::MaxNumContiguous(parsing::StringToRna("")));
  EXPECT_EQ(1, fold::MaxNumContiguous(parsing::StringToRna("A")));
  EXPECT_EQ(2, fold::MaxNumContiguous(parsing::StringToRna("AA")));
  EXPECT_EQ(2, fold::MaxNumContiguous(parsing::StringToRna("GUAAC")));
  EXPECT_EQ(1, fold::MaxNumContiguous(parsing::StringToRna("GUACA")));
  EXPECT_EQ(3, fold::MaxNumContiguous(parsing::StringToRna("GAUCCC")));
  EXPECT_EQ(3, fold::MaxNumContiguous(parsing::StringToRna("GGGAUC")));
  EXPECT_EQ(4, fold::MaxNumContiguous(parsing::StringToRna("GGGAUCAAAA")));
  EXPECT_EQ(5, fold::MaxNumContiguous(parsing::StringToRna("GGGAUUUUUCAAAA")));
}

}
}
