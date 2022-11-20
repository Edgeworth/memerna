// Copyright 2016 E.
#include "compute/mfe/mfe.h"

#include <string>

#include "common_test.h"
#include "ctx/ctx.h"
#include "ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::mfe {

class MfeAlgTest : public testing::TestWithParam<ctx::CtxCfg::DpAlg> {
 public:
  static Energy Mfe(const std::string& s) {
    return ctx::Ctx(t04, ctx::CtxCfg{.dp_alg = GetParam()}).Fold(Primary::FromSeq(s)).mfe.energy;
  }
};

TEST_P(MfeAlgTest, T04) {
  // Fast enough for brute force:
  EXPECT_EQ(E(-0.6), Mfe("CCUCCGGG"));
  EXPECT_EQ(E(-0.6), Mfe("CGGAAACGG"));
  EXPECT_EQ(E(-0.4), Mfe("UGCAAAGCAA"));
  EXPECT_EQ(E(-4.5), Mfe("GGGGAAACCCC"));
  EXPECT_EQ(E(-1.2), Mfe("CUUAUAGUUAAGG"));
  EXPECT_EQ(E(-4.0), Mfe("CCGAAGGGGCUGCGGCG"));
  EXPECT_EQ(E(-2.9), Mfe("GCCAAGGCCCCACCCGGA"));
  EXPECT_EQ(E(-4.9), Mfe("GGCCGAUGGCAGCGAUAGC"));
  EXPECT_EQ(E(-2.2), Mfe("CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::DpAlg::BRUTE) return;
  EXPECT_EQ(E(-5.1), Mfe("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  EXPECT_EQ(
      E(-13.3), Mfe("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  EXPECT_EQ(E(-5.7), Mfe("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  EXPECT_EQ(E(-12.1), Mfe("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  EXPECT_EQ(E(-7.4), Mfe("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  EXPECT_EQ(E(-3.2), Mfe("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  EXPECT_EQ(E(-12.0), Mfe("CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  EXPECT_EQ(E(-7.4), Mfe("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  EXPECT_EQ(E(-3.0), Mfe("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  EXPECT_EQ(E(-6.5), Mfe("CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  EXPECT_EQ(E(-6.0), Mfe("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  EXPECT_EQ(E(-12.2), Mfe("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  EXPECT_EQ(E(-3.9), Mfe("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  EXPECT_EQ(E(-6.7), Mfe("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  EXPECT_EQ(E(-27.6),
      Mfe("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  EXPECT_EQ(E(-5.3), Mfe("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  EXPECT_EQ(E(-15.7), Mfe("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  EXPECT_EQ(E(-4.4), Mfe("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  EXPECT_EQ(E(-2.9), Mfe("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  EXPECT_EQ(E(-2.3), Mfe("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  EXPECT_EQ(E(-8.0), Mfe("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  EXPECT_EQ(E(-20.8), Mfe("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
}

INSTANTIATE_TEST_SUITE_P(FoldAlgTest, MfeAlgTest, testing::ValuesIn(ctx::CtxCfg::DP_ALGS));

}  // namespace mrna::mfe
