// Copyright 2016 E.
#include "api/mfe.h"

#include <string>

#include "api/energy/model.h"
#include "common_test.h"
#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::mfe {

class T22MfeTest : public testing::TestWithParam<ctx::CtxCfg::DpAlg> {
 public:
  static Energy Mfe(const erg::EnergyModelPtr& em, const std::string& s) {
    return ctx::Ctx(em, ctx::CtxCfg{.dp_alg = GetParam()}).Fold(Primary::FromSeq(s)).mfe.energy;
  }
};

TEST_P(T22MfeTest, T22P2) {
  // Fast enough for brute force:
  EXPECT_EQ(E(-0.56), Mfe(t22p2, "CCUCCGGG"));
  EXPECT_EQ(E(-0.56), Mfe(t22p2, "CGGAAACGG"));
  EXPECT_EQ(E(-0.48), Mfe(t22p2, "UGCAAAGCAA"));
  EXPECT_EQ(E(-4.38), Mfe(t22p2, "GGGGAAACCCC"));
  EXPECT_EQ(E(-1.29), Mfe(t22p2, "CUUAUAGUUAAGG"));
  EXPECT_EQ(E(-3.94), Mfe(t22p2, "CCGAAGGGGCUGCGGCG"));
  EXPECT_EQ(E(-2.88), Mfe(t22p2, "GCCAAGGCCCCACCCGGA"));
  EXPECT_EQ(E(-4.90), Mfe(t22p2, "GGCCGAUGGCAGCGAUAGC"));
  EXPECT_EQ(E(-2.19), Mfe(t22p2, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::DpAlg::BRUTE) return;
  EXPECT_EQ(E(-5.25), Mfe(t22p2, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  EXPECT_EQ(E(-13.47),
      Mfe(t22p2, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  EXPECT_EQ(E(-5.69), Mfe(t22p2, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  EXPECT_EQ(E(-12.08), Mfe(t22p2, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  EXPECT_EQ(E(-7.33), Mfe(t22p2, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  EXPECT_EQ(E(-3.14), Mfe(t22p2, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  EXPECT_EQ(E(-11.90), Mfe(t22p2, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  EXPECT_EQ(E(-7.45), Mfe(t22p2, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  EXPECT_EQ(E(-2.97), Mfe(t22p2, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  EXPECT_EQ(E(-6.41), Mfe(t22p2, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  EXPECT_EQ(E(-6.07), Mfe(t22p2, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  EXPECT_EQ(
      E(-12.04), Mfe(t22p2, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  EXPECT_EQ(E(-4.10), Mfe(t22p2, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  EXPECT_EQ(E(-6.73), Mfe(t22p2, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  EXPECT_EQ(E(-27.66),
      Mfe(t22p2, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  EXPECT_EQ(E(-5.35), Mfe(t22p2, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  EXPECT_EQ(E(-15.70), Mfe(t22p2, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  EXPECT_EQ(E(-4.42), Mfe(t22p2, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  EXPECT_EQ(E(-3.06), Mfe(t22p2, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  EXPECT_EQ(E(-2.30), Mfe(t22p2, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  EXPECT_EQ(E(-8.12), Mfe(t22p2, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  EXPECT_EQ(E(-20.82), Mfe(t22p2, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
}

INSTANTIATE_TEST_SUITE_P(FoldAlgTest, T22MfeTest, testing::ValuesIn(ctx::CtxCfg::DP_ALGS));

}  // namespace mrna::mfe
