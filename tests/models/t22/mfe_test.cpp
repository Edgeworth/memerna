// Copyright 2016 Eliot Courtney.
#include "api/mfe.h"

#include <string>
#include <tuple>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/energy/model.h"
#include "api/trace/trace.h"
#include "common_test.h"
#include "gtest/gtest.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::mfe {

class T22MfeTest : public testing::TestWithParam<CtxCfg::DpAlg> {
 public:
  static std::tuple<Energy, std::string> Mfe(const erg::EnergyModelPtr& em, const std::string& s) {
    return Mfe(em, Primary::FromSeq(s));
  }

  static std::tuple<Energy, std::string> Mfe(const erg::EnergyModelPtr& em, const Primary& r) {
    auto res = Ctx(em, CtxCfg{.dp_alg = GetParam()}).Fold(r, {});
    return {res.mfe.energy, res.tb.ctd.ToString(res.tb.s)};
  }
};

#if ENERGY_PRECISION == 2

TEST_P(T22MfeTest, T22P2) {
  auto em = t22p2;

  // Fast enough for brute force:
  std::tuple<Energy, std::string> ans = {E(-0.58), "[[....]]"};
  EXPECT_EQ(ans, Mfe(em, "CCUCCGGG"));
  ans = {E(-0.53), "[[....]]3"};
  EXPECT_EQ(ans, Mfe(em, "CGGAAACGG"));
  ans = {E(-0.49), "[[[...]]]3"};
  EXPECT_EQ(ans, Mfe(em, "UGCAAAGCAA"));
  ans = {E(-4.44), "[[[[...]]]]"};
  EXPECT_EQ(ans, Mfe(em, "GGGGAAACCCC"));
  ans = {E(-1.42), "[[[[....]]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CUUAUAGUUAAGG"));
  ans = {E(-4.02), "n[....]]p[....]]."};
  EXPECT_EQ(ans, Mfe(em, "CCGAAGGGGCUGCGGCG"));
  ans = {E(-2.94), "n[....]]mp[....]]M"};
  EXPECT_EQ(ans, Mfe(em, "GCCAAGGCCCCACCCGGA"));
  ans = {E(-5.00), "mn[[...]]]Mp[....]]"};
  EXPECT_EQ(ans, Mfe(em, "GGCCGAUGGCAGCGAUAGC"));
  ans = {E(-2.08), "......[[[....]]]3...."};
  EXPECT_EQ(ans, Mfe(em, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == CtxCfg::DpAlg::BRUTE) return;
  ans = {E(-5.31), ".mn[...[[[....]]]]]Mp[[.[[[[..[[................]]..]]]]]]]."};
  EXPECT_EQ(ans, Mfe(em, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  ans = {E(-13.48), "....n[[[[...]]]]]p[[[[[............[[..[[[...]]]..]]............]]]]]].."};
  EXPECT_EQ(
      ans, Mfe(em, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  ans = {E(-5.83), "...m[[[[[[....]]]3mn[[.....]]]MP]]M"};
  EXPECT_EQ(ans, Mfe(em, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  ans = {E(-12.39), "m[[[[[3n[....]]mp[[[[.......]]]]]M]]]]]M."};
  EXPECT_EQ(ans, Mfe(em, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  ans = {E(-7.54), ".....5[[[[........[[[[[......]]]]]..]]]]"};
  EXPECT_EQ(ans, Mfe(em, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  ans = {E(-3.24), "[[[3...n[[...]]]mp[.....]]M]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  ans = {E(-12.06), ".n[[[....]]]]p[[[....]]]]..."};
  EXPECT_EQ(ans, Mfe(em, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  ans = {E(-7.57), "[[[[M..[[[......]]]3n[[....]]]mP]]]"};
  EXPECT_EQ(ans, Mfe(em, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  ans = {E(-3.75), "[[[[.[[.[[[....]]]...]].].]]]3."};
  EXPECT_EQ(ans, Mfe(em, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  ans = {E(-6.47), "[[[3n[[[...]]]]p[.....]]]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  ans = {E(-6.72), ".[[[[..[[[[....]]]].........]]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  ans = {E(-11.73), "...[[[[[M...[[.[[.[[[[...]]]].]].]]3.n[[..[[......]]..]]]mP]]]]"};
  EXPECT_EQ(ans, Mfe(em, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  ans = {E(-4.57), "[[[3..n[[[........]]]]p[[[[....]]]]]..........]]]3"};
  EXPECT_EQ(ans, Mfe(em, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  ans = {E(-5.95), "..m[[[[[[.....]]]]]]M..[[....]]3"};
  EXPECT_EQ(ans, Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  ans = {
      E(-30.43), "[[[[[[[[[[.[[[.[[.[[....[[.[[[[[[......]]].]]].]]]].]].]]]...]]]..]]]]]]]3..."};
  EXPECT_EQ(ans,
      Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  ans = {E(-5.33), "......n[[[....]]]]p[[......]]]"};
  EXPECT_EQ(ans, Mfe(em, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  ans = {E(-15.33), "[[[3....n[[[[[[[[[....]]]]...]]]]]]p[[....]]].]]]3"};
  EXPECT_EQ(ans, Mfe(em, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  ans = {E(-4.24), "[[[..........]]]3.............."};
  EXPECT_EQ(ans, Mfe(em, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  ans = {E(-3.59), "[[[3n[[[....]]]]mp[[....]]]M...]]]3"};
  EXPECT_EQ(ans, Mfe(em, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  ans = {E(-3.00), ".5[[[..[[[[....]]]]........]]]"};
  EXPECT_EQ(ans, Mfe(em, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  ans = {E(-8.44), "[[[[[3n[[....]]]p[[[...]]]]]]]]]"};
  EXPECT_EQ(ans, Mfe(em, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  ans = {E(-19.97), "[[[[[Mn[[[[[[...]]]]]]]mp[[[[[[.......]]]]]]]M.m]]]]]3..."};
  EXPECT_EQ(ans, Mfe(em, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
  ans = {E(-92.19),
      "......m[[...[[[[[[[[............[[..[[[[[[[Mp[.[[[[........]]]]..]]....n[[[[[n[[[[...[[.[[[."
      "....[[[.[[.[[[[......]]]].]].]]]......]]].]]]]]]]mp[[.[[[[....]]]]..]]]M...n[[[[[[[[[[......"
      ".......]]]]]]]]]]]mp[[[.....]]]]M...mn[...[[[[[[[.......]]]]]]]]]MP]]]]]p[[.[[[[....[[[...]]"
      "]....]]]].]]]mN]]]]]]..]]]]]]]]]]...]]M"};
  EXPECT_EQ(ans, Mfe(em, std::get<Primary>(k16sHSapiens3)));
}

#else

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(T22MfeTest);

#endif

INSTANTIATE_TEST_SUITE_P(FoldAlgTest, T22MfeTest, testing::ValuesIn(CtxCfg::DP_ALGS));

}  // namespace mrna::mfe
