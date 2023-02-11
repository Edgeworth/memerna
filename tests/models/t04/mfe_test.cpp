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

namespace mrna::md::t04 {

class T04MfeTest : public testing::TestWithParam<CtxCfg::DpAlg> {
 public:
  static std::tuple<Energy, std::string> Mfe(const erg::EnergyModelPtr& em, const std::string& s) {
    return Mfe(em, Primary::FromSeq(s));
  }

  static std::tuple<Energy, std::string> Mfe(const erg::EnergyModelPtr& em, const Primary& r) {
    auto res = Ctx(em, CtxCfg{.dp_alg = GetParam()}).Fold(r);
    return {res.mfe.energy, res.tb.ctd.ToString(res.tb.s)};
  }
};

#if ENERGY_PRECISION == 1

TEST_P(T04MfeTest, T04P1) {
  auto em = t04p1;

  // Fast enough for brute force:
  auto ans = {E(-0.6), ""};
  EXPECT_EQ(ans, Mfe(em, "CCUCCGGG"));
  ans = {E(-0.6), ""};
  EXPECT_EQ(ans, Mfe(em, "CGGAAACGG"));
  ans = {E(-0.4), ""};
  EXPECT_EQ(ans, Mfe(em, "UGCAAAGCAA"));
  ans = {E(-4.5), ""};
  EXPECT_EQ(ans, Mfe(em, "GGGGAAACCCC"));
  ans = {E(-1.2), ""};
  EXPECT_EQ(ans, Mfe(em, "CUUAUAGUUAAGG"));
  ans = {E(-4.0), ""};
  EXPECT_EQ(ans, Mfe(em, "CCGAAGGGGCUGCGGCG"));
  ans = {E(-2.9), ""};
  EXPECT_EQ(ans, Mfe(em, "GCCAAGGCCCCACCCGGA"));
  ans = {E(-4.9), ""};
  EXPECT_EQ(ans, Mfe(em, "GGCCGAUGGCAGCGAUAGC"));
  ans = {E(-2.2), ""};
  EXPECT_EQ(ans, Mfe(em, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == CtxCfg::DpAlg::BRUTE) return;
  ans = {E(-5.1), ""};
  EXPECT_EQ(ans, Mfe(em, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  ans = {E(-13.3), ""};
  EXPECT_EQ(
      ans, Mfe(em, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  ans = {E(-5.7), ""};
  EXPECT_EQ(ans, Mfe(em, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  ans = {E(-12.1), ""};
  EXPECT_EQ(ans, Mfe(em, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  ans = {E(-7.4), ""};
  EXPECT_EQ(ans, Mfe(em, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  ans = {E(-3.2), ""};
  EXPECT_EQ(ans, Mfe(em, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  ans = {E(-12.0), ""};
  EXPECT_EQ(ans, Mfe(em, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  ans = {E(-7.4), ""};
  EXPECT_EQ(ans, Mfe(em, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  ans = {E(-3.0), ""};
  EXPECT_EQ(ans, Mfe(em, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  ans = {E(-6.5), ""};
  EXPECT_EQ(ans, Mfe(em, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  ans = {E(-6.0), ""};
  EXPECT_EQ(ans, Mfe(em, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  EXPECT_EQ(ans = {E(-12.2), ""};
            EXPECT_EQ(ans, Mfe(em, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  ans = {E(-3.9), ""};
  EXPECT_EQ(ans, Mfe(em, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
ans = {E(-6.7), ""};
  EXPECT_EQ(ans, Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
ans = {E(-27.6), ""};
  EXPECT_EQ(ans,
      Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
ans = {E(-5.3), ""};
  EXPECT_EQ(ans, Mfe(em, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
ans = {E(-15.7), ""};
  EXPECT_EQ(ans, Mfe(em, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
ans = {E(-4.4), ""};
  EXPECT_EQ(ans, Mfe(em, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
ans = {E(-2.9), ""};
  EXPECT_EQ(ans, Mfe(em, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
ans = {E(-2.3), ""};
  EXPECT_EQ(ans, Mfe(em, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
ans = {E(-8.0), ""};
  EXPECT_EQ(ans, Mfe(em, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
ans = {E(-20.8), ""};
  EXPECT_EQ(ans, Mfe(em, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
ans = {E(-20.8), ""};
  EXPECT_EQ(ans, Mfe(em, std::get<Primary>(k16sHSapiens3)));
}

#elif ENERGY_PRECISION == 2

TEST_P(T04MfeTest, T04P2) {
  auto em = t04p2;

  // Fast enough for brute force:
  std::tuple<Energy, std::string> ans = {E(-0.56), "[[....]]"};
  EXPECT_EQ(ans, Mfe(em, "CCUCCGGG"));
  ans = {E(-0.56), "[[....]]3"};
  EXPECT_EQ(ans, Mfe(em, "CGGAAACGG"));
  ans = {E(-0.48), "[[[...]]]3"};
  EXPECT_EQ(ans, Mfe(em, "UGCAAAGCAA"));
  ans = {E(-4.38), "[[[[...]]]]"};
  EXPECT_EQ(ans, Mfe(em, "GGGGAAACCCC"));
  ans = {E(-1.29), "[[[[....]]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CUUAUAGUUAAGG"));
  ans = {E(-3.94), "n[....]]p[....]]."};
  EXPECT_EQ(ans, Mfe(em, "CCGAAGGGGCUGCGGCG"));
  ans = {E(-2.88), "n[....]]mp[....]]M"};
  EXPECT_EQ(ans, Mfe(em, "GCCAAGGCCCCACCCGGA"));
  ans = {E(-4.90), "mn[[...]]]Mp[....]]"};
  EXPECT_EQ(ans, Mfe(em, "GGCCGAUGGCAGCGAUAGC"));
  ans = {E(-2.19), "......[[[....]]]3...."};
  EXPECT_EQ(ans, Mfe(em, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == CtxCfg::DpAlg::BRUTE) return;
  ans = {E(-5.25), "......m[[[[...[[[..[[[...]]]...]]].]]]]M...................."};
  EXPECT_EQ(ans, Mfe(em, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  ans = {E(-13.47), "....n[[[[...]]]]]p[[[[[3...............mn[[[[[...]]]]]]Mp[....]]]]]]]].."};
  EXPECT_EQ(
      ans, Mfe(em, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  ans = {E(-5.69), "...m[[[[[[....]]]3mn[[.....]]]MP]]M"};
  EXPECT_EQ(ans, Mfe(em, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  ans = {E(-12.08), "m[[[[[3n[....]]mp[[[[.......]]]]]M]]]]]M."};
  EXPECT_EQ(ans, Mfe(em, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  ans = {E(-7.33), ".....5[[[[........[[[[[......]]]]]..]]]]"};
  EXPECT_EQ(ans, Mfe(em, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  ans = {E(-3.14), "[[[3...n[[...]]]mp[.....]]M]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  ans = {E(-11.90), ".n[[[....]]]]p[[[....]]]]..."};
  EXPECT_EQ(ans, Mfe(em, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  ans = {E(-7.45), "[[[[M..[[[......]]]3n[[....]]]mP]]]"};
  EXPECT_EQ(ans, Mfe(em, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  ans = {E(-2.97), "[[[[[[...]]]3...mn[....]]MP]]3."};
  EXPECT_EQ(ans, Mfe(em, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  ans = {E(-6.41), "[[[3n[[[...]]]]p[.....]]]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  ans = {E(-6.07), ".[[[[..[[[[....]]]].........]]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  ans = {E(-12.04), "...[[[[[M...[[.[[.[[[.....]]].]].]]3.n[[..[[......]]..]]]mP]]]]"};
  EXPECT_EQ(ans, Mfe(em, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  ans = {E(-4.10), "[[[3..n[[[........]]]]p[[[[....]]]]]..........]]]3"};
  EXPECT_EQ(ans, Mfe(em, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  ans = {E(-6.73), "..m[[[[[[.....]]]]]]M..[[....]]3"};
  EXPECT_EQ(ans, Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  ans = {
      E(-27.66), "[[[[[[[[[[[.[[...[[[[....]]]]..]].]]]3mn[[[..[[[[....]]]]...]]]]MP]]]]]]]3..."};
  EXPECT_EQ(ans,
      Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  ans = {E(-5.35), "......n[[[....]]]]p[[......]]]"};
  EXPECT_EQ(ans, Mfe(em, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  ans = {E(-15.70), "[[[3....n[[[[[[[[[....]]]]...]]]]]]p[[....]]].]]]3"};
  EXPECT_EQ(ans, Mfe(em, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  ans = {E(-4.42), ".[[[[3.n[....]]mp[....]]M..]]]]"};
  EXPECT_EQ(ans, Mfe(em, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  ans = {E(-3.06), "[[[3n[[[....]]]]mp[[.....]]]M..]]]3"};
  EXPECT_EQ(ans, Mfe(em, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  ans = {E(-2.30), ".5[[[[[[...]]]3mn[[....]]]MP]]"};
  EXPECT_EQ(ans, Mfe(em, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  ans = {E(-8.12), "m[[[[[[[[....]]]3n[[...]]]P]]]]M"};
  EXPECT_EQ(ans, Mfe(em, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  ans = {E(-20.82), "[[[[[Mn[[[[[[...]]]]]]]mp[[[[[[.......]]]]]]]M.m]]]]]3..."};
  EXPECT_EQ(ans, Mfe(em, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
  ans = {E(-89.99),
      "......m[[...[[[[[[[[............[[..[[[[[[[M..m[[[[........]]]]M.......[[[[[[n[[[[...[[.[[[."
      "....[[[.[[.[[[[......]]]].]].]]]......]]].]]]]]]]mp[[.[[[[....]]]]..]]]M...n[[[[[[[[[[......"
      ".......]]]]]]]]]]]mp[[[.....]]]]M...mn[...[[[[[[[.......]]]]]]]]]MP]]]]]3...........n[[[[[.."
      ".......]]]]]]mP]]]]]]..]]]]]]]]]]...]]M"};
  EXPECT_EQ(ans, Mfe(em, std::get<Primary>(k16sHSapiens3)));
}

TEST_P(T04MfeTest, T12P2) {
  auto em = t12p2;

  // Fast enough for brute force:
  std::tuple<Energy, std::string> ans = {E(-0.56), "[[....]]"};
  EXPECT_EQ(ans, Mfe(em, "CCUCCGGG"));
  ans = {E(-0.56), "[[....]]3"};
  EXPECT_EQ(ans, Mfe(em, "CGGAAACGG"));
  ans = {E(-0.48), "[[[...]]]3"};
  EXPECT_EQ(ans, Mfe(em, "UGCAAAGCAA"));
  ans = {E(-4.38), "[[[[...]]]]"};
  EXPECT_EQ(ans, Mfe(em, "GGGGAAACCCC"));
  ans = {E(-1.29), "[[[[....]]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CUUAUAGUUAAGG"));
  ans = {E(-3.94), "n[....]]p[....]]."};
  EXPECT_EQ(ans, Mfe(em, "CCGAAGGGGCUGCGGCG"));
  ans = {E(-2.88), "n[....]]mp[....]]M"};
  EXPECT_EQ(ans, Mfe(em, "GCCAAGGCCCCACCCGGA"));
  ans = {E(-4.90), "mn[[...]]]Mp[....]]"};
  EXPECT_EQ(ans, Mfe(em, "GGCCGAUGGCAGCGAUAGC"));
  ans = {E(-2.19), "......[[[....]]]3...."};
  EXPECT_EQ(ans, Mfe(em, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (GetParam() == CtxCfg::DpAlg::BRUTE) return;
  ans = {E(-5.34), "......m[[[[...[[[..[[[...]]]...]]].]]]]M...................."};
  EXPECT_EQ(ans, Mfe(em, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  ans = {E(-13.09), "....n[[[[...]]]]]p[[[[[............[[..[[[...]]]..]]............]]]]]].."};
  EXPECT_EQ(
      ans, Mfe(em, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  ans = {E(-5.69), "...m[[[[[[....]]]3mn[[.....]]]MP]]M"};
  EXPECT_EQ(ans, Mfe(em, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  ans = {E(-12.08), "m[[[[[3n[....]]mp[[[[.......]]]]]M]]]]]M."};
  EXPECT_EQ(ans, Mfe(em, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  ans = {E(-7.33), ".....5[[[[........[[[[[......]]]]]..]]]]"};
  EXPECT_EQ(ans, Mfe(em, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  ans = {E(-2.83), "[[[3...n[[...]]]mp[.....]]M]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  ans = {E(-11.83), ".n[[[....]]]]p[[[....]]]]..."};
  EXPECT_EQ(ans, Mfe(em, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  ans = {E(-7.54), "[[[[M..[[[......]]]3n[[....]]]mP]]]"};
  EXPECT_EQ(ans, Mfe(em, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  ans = {E(-3.18), "........m[[[.....[[....]]..]]]M"};
  EXPECT_EQ(ans, Mfe(em, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  ans = {E(-6.41), "[[[3n[[[...]]]]p[.....]]]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  ans = {E(-5.93), "m[[[...]]]M.....[[[[[....]]..]]]3"};
  EXPECT_EQ(ans, Mfe(em, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  ans = {E(-12.10), "...[[[[[M...[[.[[.[[[[...]]]].]].]]3.n[[..[[......]]..]]]mP]]]]"};
  EXPECT_EQ(ans, Mfe(em, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  ans = {E(-4.10), "[[[3..n[[[........]]]]p[[[[....]]]]]..........]]]3"};
  EXPECT_EQ(ans, Mfe(em, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  ans = {E(-5.98), "..m[[[[[[.....]]]]]]M..[[....]]3"};
  EXPECT_EQ(ans, Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  ans = {
      E(-28.39), "[[[[[[[[[[[.[[...[[[[....]]]]..]].]]]3mn[[[..[[[[....]]]]...]]]]MP]]]]]]]3..."};
  EXPECT_EQ(ans,
      Mfe(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  ans = {E(-5.44), "......n[[[....]]]]p[[......]]]"};
  EXPECT_EQ(ans, Mfe(em, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  ans = {E(-15.11), "[[[3....n[[[[[[[[[....]]]]...]]]]]]p[[....]]].]]]3"};
  EXPECT_EQ(ans, Mfe(em, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  ans = {E(-4.18), "[[[..........]]]3.............."};
  EXPECT_EQ(ans, Mfe(em, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  ans = {E(-3.35), "[[[3n[[[....]]]]mp[[.....]]]M..]]]3"};
  EXPECT_EQ(ans, Mfe(em, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  ans = {E(-2.39), ".5[[[[[[...]]]3mn[[....]]]MP]]"};
  EXPECT_EQ(ans, Mfe(em, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  ans = {E(-8.53), "[[[[[[[[[....]]]3n[[...]]]P]]]]]"};
  EXPECT_EQ(ans, Mfe(em, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  ans = {E(-20.14), "[[[[[Mn[[[[[[...]]]]]]]mp[[[[[[.......]]]]]]]M.m]]]]]3..."};
  EXPECT_EQ(ans, Mfe(em, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
  ans = {E(-89.54),
      "......m[[...[[[[[[[[............[[..[[[[[[[M..m[[[[........]]]]M.......[[[[[[n[[[[...[[.[[[."
      "....[[[.[[.[[[[......]]]].]].]]]......]]].]]]]]]]mp[[.[[[[....]]]]..]]]M...n[[[[[[[[[[......"
      ".......]]]]]]]]]]]mp[[[.....]]]]M...mn[...[[[[[[[.......]]]]]]]]]MP]]]]]3...........n[[[[[.."
      ".......]]]]]]mP]]]]]]..]]]]]]]]]]...]]M"};
  EXPECT_EQ(ans, Mfe(em, std::get<Primary>(k16sHSapiens3)));
}

// NEWMODEL: Add tests here.

#endif

INSTANTIATE_TEST_SUITE_P(FoldAlgTest, T04MfeTest, testing::ValuesIn(CtxCfg::DP_ALGS));

}  // namespace mrna::md::t04
