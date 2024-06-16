// Copyright 2024 Eliot Courtney.
#include <string>
#include <tuple>

#include "api/ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"
#include "tests/init.h"
#include "tests/util.h"

namespace mrna {

class MfeTestT04 : public testing::TestWithParam<std::tuple<int, CtxCfg::MfeAlg>> {
 public:
  static std::tuple<Energy, std::string> Mfe(const BackendModelPtr& m, const std::string& s) {
    return GetMfe(m, std::get<1>(GetParam()), s);
  }

  static std::tuple<Energy, std::string> Mfe(const BackendModelPtr& m, const Primary& r) {
    return GetMfe(m, std::get<1>(GetParam()), r);
  }
};

#if ENERGY_PRECISION == 1

TEST_P(MfeTestT04, T04P1) {
  auto [i, alg] = GetParam();
  auto m = t04_ms[i];
  if (!Contains(CtxCfg::MfeAlgsForBackend(m), alg)) return;

  // Fast enough for brute force:
  std::tuple<Energy, std::string> ans = {E(-0.6), "[[....]]"};
  EXPECT_EQ(ans, Mfe(m, "CCUCCGGG"));
  ans = {E(-0.6), "[[....]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGGAAACGG"));
  ans = {E(-0.4), "[[[...]]]3"};
  EXPECT_EQ(ans, Mfe(m, "UGCAAAGCAA"));
  ans = {E(-4.5), "[[[[...]]]]"};
  EXPECT_EQ(ans, Mfe(m, "GGGGAAACCCC"));
  ans = {E(-1.2), "[[[[....]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CUUAUAGUUAAGG"));
  ans = {E(-4.0), "n[....]]p[....]]."};
  EXPECT_EQ(ans, Mfe(m, "CCGAAGGGGCUGCGGCG"));
  ans = {E(-2.9), "n[....]]mp[....]]M"};
  EXPECT_EQ(ans, Mfe(m, "GCCAAGGCCCCACCCGGA"));
  ans = {E(-4.9), "mn[[...]]]Mp[....]]"};
  EXPECT_EQ(ans, Mfe(m, "GGCCGAUGGCAGCGAUAGC"));
  ans = {E(-2.2), "......[[[....]]]3...."};
  EXPECT_EQ(ans, Mfe(m, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (alg == CtxCfg::MfeAlg::BRUTE) return;
  ans = {E(-5.1), "......m[[[[...[[[..[[[...]]]...]]].]]]]M...................."};
  EXPECT_EQ(ans, Mfe(m, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  ans = {E(-13.3), ".....n[[[...]]]]mp[[[[[3...............mn[[[[[...]]]]]]Mp[....]]]]]]]]M."};
  EXPECT_EQ(
      ans, Mfe(m, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  ans = {E(-5.7), "...m[[[[[[....]]]3mn[[.....]]]MP]]M"};
  EXPECT_EQ(ans, Mfe(m, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  ans = {E(-12.1), "m[[[[[3n[....]]mp[[[[.......]]]]]M]]]]]M."};
  EXPECT_EQ(ans, Mfe(m, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  ans = {E(-7.4), ".....5[[[[........[[[[[......]]]]]..]]]]"};
  EXPECT_EQ(ans, Mfe(m, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  ans = {E(-3.2), "[[[3...n[[...]]]mp[.....]]M]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  ans = {E(-12.0), ".n[[[....]]]]p[[[....]]]]..."};
  EXPECT_EQ(ans, Mfe(m, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  ans = {E(-7.4), "[[[[M..[[[......]]]3n[[....]]]mP]]]"};
  EXPECT_EQ(ans, Mfe(m, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  ans = {E(-3.0), "[[[[[[...]]]3...mn[....]]MP]]3."};
  EXPECT_EQ(ans, Mfe(m, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  ans = {E(-6.5), "[[[3n[[[...]]]]p[.....]]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  ans = {E(-6.0), ".[[[[..[[[[....]]]].........]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  ans = {E(-12.2), "...[[[[[M...[[.[[.[[[.....]]].]].]]3.n[[..[[......]]..]]]mP]]]]"};
  EXPECT_EQ(ans, Mfe(m, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  ans = {E(-3.9), "[[[3..n[[[........]]]]p[[[[....]]]]]..........]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  ans = {E(-6.7), "..m[[[[[[.....]]]]]]M..[[....]]3"};
  EXPECT_EQ(ans, Mfe(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  ans = {E(-27.6), "[[[[[[[[[[[.[[...[[[[....]]]]..]].]]]3mn[[[..[[[[....]]]]...]]]]MP]]]]]]]3..."};
  EXPECT_EQ(
      ans, Mfe(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  ans = {E(-5.3), "......n[[[....]]]]p[[......]]]"};
  EXPECT_EQ(ans, Mfe(m, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  ans = {E(-15.7), "[[[3....n[[[[[[[[[....]]]]...]]]]]]p[[....]]].]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  ans = {E(-4.4), ".[[[[3.n[....]]mp[....]]M..]]]]"};
  EXPECT_EQ(ans, Mfe(m, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  ans = {E(-2.9), "[[[3n[[[....]]]]mp[[.....]]]M..]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  ans = {E(-2.3), ".5[[[[[[...]]]3mn[[....]]]MP]]"};
  EXPECT_EQ(ans, Mfe(m, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));

  // Has multiple possible MFE structures.
  EXPECT_EQ(E(-8.0), std::get<0>(Mfe(m, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG")));

  ans = {E(-20.8), ".[[[[Mn[[[[[[...]]]]]]]mp[[[[[[.......]]]]]]]M.m]]]]3...."};
  EXPECT_EQ(ans, Mfe(m, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));

  // Has multiple possible MFE structures.
  EXPECT_EQ(E(-89.6), std::get<0>(Mfe(m, std::get<Primary>(k16sHSapiens3))));
}

#elif ENERGY_PRECISION == 2

TEST_P(MfeTestT04, T04P2) {
  auto [i, alg] = GetParam();
  auto m = t04_ms[i];
  if (!Contains(CtxCfg::MfeAlgsForBackend(m), alg)) return;

  // Fast enough for brute force:
  std::tuple<Energy, std::string> ans = {E(-0.56), "[[....]]"};
  EXPECT_EQ(ans, Mfe(m, "CCUCCGGG"));
  ans = {E(-0.56), "[[....]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGGAAACGG"));
  ans = {E(-0.43), "[[[...]]]3"};
  EXPECT_EQ(ans, Mfe(m, "UGCAAAGCAA"));
  ans = {E(-4.38), "[[[[...]]]]"};
  EXPECT_EQ(ans, Mfe(m, "GGGGAAACCCC"));
  ans = {E(-1.24), "[[[[....]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CUUAUAGUUAAGG"));
  ans = {E(-3.94), "n[....]]p[....]]."};
  EXPECT_EQ(ans, Mfe(m, "CCGAAGGGGCUGCGGCG"));
  ans = {E(-2.88), "n[....]]mp[....]]M"};
  EXPECT_EQ(ans, Mfe(m, "GCCAAGGCCCCACCCGGA"));
  ans = {E(-4.90), "mn[[...]]]Mp[....]]"};
  EXPECT_EQ(ans, Mfe(m, "GGCCGAUGGCAGCGAUAGC"));
  ans = {E(-2.19), "......[[[....]]]3...."};
  EXPECT_EQ(ans, Mfe(m, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (alg == CtxCfg::MfeAlg::BRUTE) return;
  ans = {E(-5.05), "......m[[[[...[[[..[[[...]]]...]]].]]]]M...................."};
  EXPECT_EQ(ans, Mfe(m, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  ans = {E(-13.32), "....n[[[[...]]]]]p[[[[[3...............mn[[[[[...]]]]]]Mp[....]]]]]]]].."};
  EXPECT_EQ(
      ans, Mfe(m, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  ans = {E(-5.59), "...m[[[[[[....]]]3mn[[.....]]]MP]]M"};
  EXPECT_EQ(ans, Mfe(m, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  ans = {E(-12.03), "m[[[[[3n[....]]mp[[[[.......]]]]]M]]]]]M."};
  EXPECT_EQ(ans, Mfe(m, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  ans = {E(-7.28), ".....5[[[[........[[[[[......]]]]]..]]]]"};
  EXPECT_EQ(ans, Mfe(m, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  ans = {E(-3.09), "[[[3...n[[...]]]mp[.....]]M]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  ans = {E(-11.90), ".n[[[....]]]]p[[[....]]]]..."};
  EXPECT_EQ(ans, Mfe(m, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  ans = {E(-7.35), "[[[[M..[[[......]]]3n[[....]]]mP]]]"};
  EXPECT_EQ(ans, Mfe(m, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  ans = {E(-2.87), "[[[[[[...]]]3...mn[....]]MP]]3."};
  EXPECT_EQ(ans, Mfe(m, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  ans = {E(-6.36), "[[[3n[[[...]]]]p[.....]]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  ans = {E(-5.97), ".[[[[..[[[[....]]]].........]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  ans = {E(-11.99), "...[[[[[M...[[.[[.[[[.....]]].]].]]3.n[[..[[......]]..]]]mP]]]]"};
  EXPECT_EQ(ans, Mfe(m, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  ans = {E(-3.95), "[[[3..n[[[........]]]]p[[[[....]]]]]..........]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  ans = {E(-6.68), "..m[[[[[[.....]]]]]]M..[[....]]3"};
  EXPECT_EQ(ans, Mfe(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  ans = {
      E(-27.41), "[[[[[[[[[[[.[[...[[[[....]]]]..]].]]]3mn[[[..[[[[....]]]]...]]]]MP]]]]]]]3..."};
  EXPECT_EQ(
      ans, Mfe(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  ans = {E(-5.25), "......n[[[....]]]]p[[......]]]"};
  EXPECT_EQ(ans, Mfe(m, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  ans = {E(-15.65), "[[[3....n[[[[[[[[[....]]]]...]]]]]]p[[....]]].]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  ans = {E(-4.42), ".[[[[3.n[....]]mp[....]]M..]]]]"};
  EXPECT_EQ(ans, Mfe(m, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  ans = {E(-2.91), "[[[3n[[[....]]]]mp[[.....]]]M..]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  ans = {E(-2.15), ".5[[[[[[...]]]3mn[[....]]]MP]]"};
  EXPECT_EQ(ans, Mfe(m, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));

  // Has multiple possible MFE structures.
  EXPECT_EQ(E(-8.07), std::get<0>(Mfe(m, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG")));

  ans = {E(-20.62), "[[[[[Mn[[[[[[...]]]]]]]mp[[[[[[.......]]]]]]]M.m]]]]]3..."};
  EXPECT_EQ(ans, Mfe(m, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
  ans = {E(-89.29),
      "......m[[...[[[[[[[[............[[..[[[[[[[M..m[[[[........]]]]M.......[[[[[[n[[[[...[[.[[[."
      "....[[[.[[.[[[[......]]]].]].]]]......]]].]]]]]]]mp[[.[[[[....]]]]..]]]M...n[[[[[[[[[[......"
      ".......]]]]]]]]]]]mp[[[.....]]]]M...mn[...[[[[[[[.......]]]]]]]]]MP]]]]]3...........n[[[[[.."
      ".......]]]]]]mP]]]]]]..]]]]]]]]]]...]]M"};
  EXPECT_EQ(ans, Mfe(m, std::get<Primary>(k16sHSapiens3)));
}

#endif

INSTANTIATE_TEST_SUITE_P(MfeTest, MfeTestT04,
    testing::Combine(
        testing::Range(0, NUM_T04_MODELS), testing::ValuesIn(EnumValues<CtxCfg::MfeAlg>())));

}  // namespace mrna
