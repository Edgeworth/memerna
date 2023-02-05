// Copyright 2016 E.
#include <string>

#include "common_test.h"
#include "compute/energy/model.h"
#include "compute/mfe/mfe.h"
#include "ctx/ctx.h"
#include "ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::mfe {

class T04MfeTest : public testing::TestWithParam<ctx::CtxCfg::DpAlg> {
 public:
  static Energy Mfe(const erg::EnergyModelPtr& em, const std::string& s) {
    return ctx::Ctx(em, ctx::CtxCfg{.dp_alg = GetParam()}).Fold(Primary::FromSeq(s)).mfe.energy;
  }
};

#if ENERGY_PRECISION == 1

TEST_P(T04MfeTest, T04P1) {
  // Fast enough for brute force:
  EXPECT_EQ(E(-0.6), Mfe(t04p1, "CCUCCGGG"));
  EXPECT_EQ(E(-0.6), Mfe(t04p1, "CGGAAACGG"));
  EXPECT_EQ(E(-0.4), Mfe(t04p1, "UGCAAAGCAA"));
  EXPECT_EQ(E(-4.5), Mfe(t04p1, "GGGGAAACCCC"));
  EXPECT_EQ(E(-1.2), Mfe(t04p1, "CUUAUAGUUAAGG"));
  EXPECT_EQ(E(-4.0), Mfe(t04p1, "CCGAAGGGGCUGCGGCG"));
  EXPECT_EQ(E(-2.9), Mfe(t04p1, "GCCAAGGCCCCACCCGGA"));
  EXPECT_EQ(E(-4.9), Mfe(t04p1, "GGCCGAUGGCAGCGAUAGC"));
  EXPECT_EQ(E(-2.2), Mfe(t04p1, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::DpAlg::BRUTE) return;
  EXPECT_EQ(E(-5.1), Mfe(t04p1, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  EXPECT_EQ(E(-13.3),
      Mfe(t04p1, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  EXPECT_EQ(E(-5.7), Mfe(t04p1, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  EXPECT_EQ(E(-12.1), Mfe(t04p1, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  EXPECT_EQ(E(-7.4), Mfe(t04p1, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  EXPECT_EQ(E(-3.2), Mfe(t04p1, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  EXPECT_EQ(E(-12.0), Mfe(t04p1, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  EXPECT_EQ(E(-7.4), Mfe(t04p1, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  EXPECT_EQ(E(-3.0), Mfe(t04p1, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  EXPECT_EQ(E(-6.5), Mfe(t04p1, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  EXPECT_EQ(E(-6.0), Mfe(t04p1, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  EXPECT_EQ(
      E(-12.2), Mfe(t04p1, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  EXPECT_EQ(E(-3.9), Mfe(t04p1, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  EXPECT_EQ(E(-6.7), Mfe(t04p1, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  EXPECT_EQ(E(-27.6),
      Mfe(t04p1, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  EXPECT_EQ(E(-5.3), Mfe(t04p1, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  EXPECT_EQ(E(-15.7), Mfe(t04p1, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  EXPECT_EQ(E(-4.4), Mfe(t04p1, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  EXPECT_EQ(E(-2.9), Mfe(t04p1, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  EXPECT_EQ(E(-2.3), Mfe(t04p1, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  EXPECT_EQ(E(-8.0), Mfe(t04p1, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  EXPECT_EQ(E(-20.8), Mfe(t04p1, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
}

#elif ENERGY_PRECISION == 2

TEST_P(T04MfeTest, T04P2) {
  // Fast enough for brute force:
  EXPECT_EQ(E(-0.56), Mfe(t04p2, "CCUCCGGG"));
  EXPECT_EQ(E(-0.56), Mfe(t04p2, "CGGAAACGG"));
  EXPECT_EQ(E(-0.48), Mfe(t04p2, "UGCAAAGCAA"));
  EXPECT_EQ(E(-4.38), Mfe(t04p2, "GGGGAAACCCC"));
  EXPECT_EQ(E(-1.29), Mfe(t04p2, "CUUAUAGUUAAGG"));
  EXPECT_EQ(E(-3.94), Mfe(t04p2, "CCGAAGGGGCUGCGGCG"));
  EXPECT_EQ(E(-2.88), Mfe(t04p2, "GCCAAGGCCCCACCCGGA"));
  EXPECT_EQ(E(-4.90), Mfe(t04p2, "GGCCGAUGGCAGCGAUAGC"));
  EXPECT_EQ(E(-2.19), Mfe(t04p2, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::DpAlg::BRUTE) return;
  EXPECT_EQ(E(-5.25), Mfe(t04p2, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  EXPECT_EQ(E(-13.47),
      Mfe(t04p2, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  EXPECT_EQ(E(-5.69), Mfe(t04p2, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  EXPECT_EQ(E(-12.08), Mfe(t04p2, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  EXPECT_EQ(E(-7.33), Mfe(t04p2, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  EXPECT_EQ(E(-3.14), Mfe(t04p2, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  EXPECT_EQ(E(-11.90), Mfe(t04p2, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  EXPECT_EQ(E(-7.45), Mfe(t04p2, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  EXPECT_EQ(E(-2.97), Mfe(t04p2, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  EXPECT_EQ(E(-6.41), Mfe(t04p2, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  EXPECT_EQ(E(-6.07), Mfe(t04p2, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  EXPECT_EQ(
      E(-12.04), Mfe(t04p2, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  EXPECT_EQ(E(-4.10), Mfe(t04p2, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  EXPECT_EQ(E(-6.73), Mfe(t04p2, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  EXPECT_EQ(E(-27.66),
      Mfe(t04p2, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  EXPECT_EQ(E(-5.35), Mfe(t04p2, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  EXPECT_EQ(E(-15.70), Mfe(t04p2, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  EXPECT_EQ(E(-4.42), Mfe(t04p2, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  EXPECT_EQ(E(-3.06), Mfe(t04p2, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  EXPECT_EQ(E(-2.30), Mfe(t04p2, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  EXPECT_EQ(E(-8.12), Mfe(t04p2, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  EXPECT_EQ(E(-20.82), Mfe(t04p2, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
}

// NEWMODEL: Add tests here.

#endif

INSTANTIATE_TEST_SUITE_P(FoldAlgTest, T04MfeTest, testing::ValuesIn(ctx::CtxCfg::DP_ALGS));

}  // namespace mrna::mfe
