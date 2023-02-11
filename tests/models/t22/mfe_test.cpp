// Copyright 2016 Eliot Courtney.
#include "api/mfe.h"

#include <string>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/energy/model.h"
#include "common_test.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::mfe {

class T22MfeTest : public testing::TestWithParam<CtxCfg::DpAlg> {
 public:
  static std::tuple<Energy, std::string> Mfe(const erg::EnergyModelPtr& em, const std::string& s) {
    return Mfe(em, Primary::FromSeq(s));
  }

  static std::tuple<Energy, std::string> Mfe(const erg::EnergyModelPtr& em, const Primary& r) {
    auto res = Ctx(em, CtxCfg{.dp_alg = GetParam()}).Fold(r);
    return {res.mfe.energy, res.tb.ctd.ToString(res.tb.s)};
  }
};

#if ENERGY_PRECISION == 2

TEST_P(T22MfeTest, T22P2) {
  auto em = t22p2;

  // Fast enough for brute force:
  std::tuple<Energy, std::string> ans = {E(-0.58), ""};
  EXPECT_EQ(ans, Mfe(em, "CCUCCGGG"));
  ans = {E(-0.53), ""};
  EXPECT_EQ(ans, Mfe(em, "CGGAAACGG"));
  ans = {E(-0.49), ""};
  EXPECT_EQ(ans, Mfe(em, "UGCAAAGCAA"));
  ans = {E(-4.44), ""};
  EXPECT_EQ(ans, Mfe(em, "GGGGAAACCCC"));
  ans = {E(-1.42), ""};
  EXPECT_EQ(ans, Mfe(em, "CUUAUAGUUAAGG"));
  ans = {E(-4.02), ""};
  EXPECT_EQ(ans, Mfe(em, "CCGAAGGGGCUGCGGCG"));
  ans = {E(-2.94), ""};
  EXPECT_EQ(ans, Mfe(em, "GCCAAGGCCCCACCCGGA"));
  ans = {E(-5.00), ""};
  EXPECT_EQ(ans, Mfe(em, "GGCCGAUGGCAGCGAUAGC"));
  ans = {E(-2.08), ""};
  EXPECT_EQ(ans, Mfe(em, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == CtxCfg::DpAlg::BRUTE) return;
  ans = {E(-5.26), ""};
  EXPECT_EQ(ans, Mfe(em, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  ans = {E(-13.43), ""};
  EXPECT_EQ(
      ans, Mfe(em, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  ans = {E(-5.83), ""};
  EXPECT_EQ(ans, Mfe(em, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  ans = {E(-12.39), ""};
  EXPECT_EQ(ans, Mfe(em, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  ans = {E(-7.54), ""};
  EXPECT_EQ(ans, Mfe(em, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  ans = {E(-3.24), ""};
  EXPECT_EQ(ans, Mfe(em, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  ans = {E(-12.06), ""};
  EXPECT_EQ(ans, Mfe(em, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  ans = {E(-7.57), ""};
  EXPECT_EQ(ans, Mfe(em, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  ans = {E(-3.70), ""};
  EXPECT_EQ(ans, Mfe(em, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  ans = {E(-6.47), ""};
  EXPECT_EQ(ans, Mfe(em, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  ans = {E(-6.72), ""};
  EXPECT_EQ(ans, Mfe(em, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  ans = {E(-11.68), ""};
  EXPECT_EQ(ans, Mfe(em, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  ans = {E(-4.57), ""};
  EXPECT_EQ(ans, Mfe(em, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  ans = {E(-5.95), ""};
  EXPECT_EQ(ans, Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  ans = {E(-30.18), ""};
  EXPECT_EQ(ans,
      Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  ans = {E(-5.33), ""};
  EXPECT_EQ(ans, Mfe(em, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  ans = {E(-15.33), ""};
  EXPECT_EQ(ans, Mfe(em, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  ans = {E(-4.24), ""};
  EXPECT_EQ(ans, Mfe(em, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  ans = {E(-3.59), ""};
  EXPECT_EQ(ans, Mfe(em, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  ans = {E(-3.00), ""};
  EXPECT_EQ(ans, Mfe(em, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  ans = {E(-8.44), ""};
  EXPECT_EQ(ans, Mfe(em, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  ans = {E(-19.97), ""};
  EXPECT_EQ(ans, Mfe(em, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));

  ans = {E(-20.8), ""};
  EXPECT_EQ(ans, Mfe(em, std::get<Primary>(k16sHSapiens3)));
}

#else

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(T22MfeTest);

#endif

INSTANTIATE_TEST_SUITE_P(FoldAlgTest, T22MfeTest, testing::ValuesIn(CtxCfg::DP_ALGS));

}  // namespace mrna::mfe
