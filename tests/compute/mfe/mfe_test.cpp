// Copyright 2016 E.
#include "compute/mfe/mfe.h"

#include <string>

#include "common_test.h"
#include "ctx/config.h"
#include "ctx/ctx.h"
#include "gtest/gtest.h"
#include "model/model.h"
#include "model/primary.h"

namespace mrna::mfe {

class MfeAlgTest : public testing::TestWithParam<ctx::CtxCfg::DpAlg> {
 public:
  Energy Mfe(const std::string& s) {
    return ctx::Ctx(t04, ctx::CtxCfg{.dp_alg = GetParam()}).Fold(Primary::FromString(s)).mfe.energy;
  }
};

TEST_P(MfeAlgTest, T04) {
  // Fast enough for brute force:
  EXPECT_EQ(-6, Mfe("CCUCCGGG"));
  EXPECT_EQ(-6, Mfe("CGGAAACGG"));
  EXPECT_EQ(-4, Mfe("UGCAAAGCAA"));
  EXPECT_EQ(-45, Mfe("GGGGAAACCCC"));
  EXPECT_EQ(-12, Mfe("CUUAUAGUUAAGG"));
  EXPECT_EQ(-40, Mfe("CCGAAGGGGCUGCGGCG"));
  EXPECT_EQ(-29, Mfe("GCCAAGGCCCCACCCGGA"));
  EXPECT_EQ(-49, Mfe("GGCCGAUGGCAGCGAUAGC"));
  EXPECT_EQ(-22, Mfe("CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::DpAlg::BRUTE) return;
  EXPECT_EQ(-51, Mfe("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  EXPECT_EQ(-133, Mfe("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  EXPECT_EQ(-57, Mfe("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  EXPECT_EQ(-121, Mfe("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  EXPECT_EQ(-74, Mfe("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  EXPECT_EQ(-32, Mfe("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  EXPECT_EQ(-120, Mfe("CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  EXPECT_EQ(-74, Mfe("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  EXPECT_EQ(-30, Mfe("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  EXPECT_EQ(-65, Mfe("CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  EXPECT_EQ(-60, Mfe("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  EXPECT_EQ(-122, Mfe("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  EXPECT_EQ(-39, Mfe("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  EXPECT_EQ(-67, Mfe("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  EXPECT_EQ(
      -276, Mfe("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  EXPECT_EQ(-53, Mfe("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  EXPECT_EQ(-157, Mfe("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  EXPECT_EQ(-44, Mfe("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  EXPECT_EQ(-29, Mfe("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  EXPECT_EQ(-23, Mfe("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  EXPECT_EQ(-80, Mfe("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  EXPECT_EQ(-208, Mfe("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
}

INSTANTIATE_TEST_SUITE_P(FoldAlgTest, MfeAlgTest, testing::ValuesIn(ctx::CtxCfg::DP_ALGS));

}  // namespace mrna::mfe
