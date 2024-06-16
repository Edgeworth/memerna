// Copyright 2024 Eliot Courtney.
#include <string>
#include <tuple>

#include "api/ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"
#include "tests/init.h"
#include "tests/util.h"
#include "util/string.h"

namespace mrna {

class MfeTestT12 : public testing::TestWithParam<std::tuple<int, CtxCfg::MfeAlg>> {
 public:
  static std::tuple<Energy, std::string> Mfe(const BackendModelPtr& m, const std::string& s) {
    return GetMfe(m, std::get<1>(GetParam()), s);
  }

  static std::tuple<Energy, std::string> Mfe(const BackendModelPtr& m, const Primary& r) {
    return GetMfe(m, std::get<1>(GetParam()), r);
  }
};

#if ENERGY_PRECISION == 1

// TODO(2): Add tests for all GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST.
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(MfeTestT12);

#elif ENERGY_PRECISION == 2

TEST_P(MfeTestT12, T12P2) {
  auto [i, alg] = GetParam();
  auto m = t12_ms[i];
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
  ans = {E(-5.19), "......m[[[[...[[[..[[[...]]]...]]].]]]]M...................."};
  EXPECT_EQ(ans, Mfe(m, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  ans = {E(-12.94), "....n[[[[...]]]]]p[[[[[............[[..[[[...]]]..]]............]]]]]].."};
  EXPECT_EQ(
      ans, Mfe(m, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  ans = {E(-5.59), "...m[[[[[[....]]]3mn[[.....]]]MP]]M"};
  EXPECT_EQ(ans, Mfe(m, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  ans = {E(-12.03), "m[[[[[3n[....]]mp[[[[.......]]]]]M]]]]]M."};
  EXPECT_EQ(ans, Mfe(m, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  ans = {E(-7.28), ".....5[[[[........[[[[[......]]]]]..]]]]"};
  EXPECT_EQ(ans, Mfe(m, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  ans = {E(-2.83), "[[[3...n[[...]]]mp[.....]]M]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  ans = {E(-11.83), ".n[[[....]]]]p[[[....]]]]..."};
  EXPECT_EQ(ans, Mfe(m, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  ans = {E(-7.49), "[[[[M..[[[......]]]3n[[....]]]mP]]]"};
  EXPECT_EQ(ans, Mfe(m, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  ans = {E(-3.18), "........m[[[.....[[....]]..]]]M"};
  EXPECT_EQ(ans, Mfe(m, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  ans = {E(-6.36), "[[[3n[[[...]]]]p[.....]]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  ans = {E(-5.93), "m[[[...]]]M.....[[[[[....]]..]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  ans = {E(-12.05), "...[[[[[M...[[.[[.[[[[...]]]].]].]]3.n[[..[[......]]..]]]mP]]]]"};
  EXPECT_EQ(ans, Mfe(m, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  ans = {E(-3.95), "[[[3..n[[[........]]]]p[[[[....]]]]]..........]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  ans = {E(-5.93), "..m[[[[[[.....]]]]]]M..[[....]]3"};
  EXPECT_EQ(ans, Mfe(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  ans = {
      E(-28.29), "[[[[[[[[[[[.[[...[[[[....]]]]..]].]]]3mn[[[..[[[[....]]]]...]]]]MP]]]]]]]3..."};
  EXPECT_EQ(
      ans, Mfe(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  ans = {E(-5.39), "......n[[[....]]]]p[[......]]]"};
  EXPECT_EQ(ans, Mfe(m, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  ans = {E(-15.06), "[[[3....n[[[[[[[[[....]]]]...]]]]]]p[[....]]].]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  ans = {E(-4.18), "[[[..........]]]3.............."};
  EXPECT_EQ(ans, Mfe(m, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  ans = {E(-3.25), "[[[3n[[[....]]]]mp[[.....]]]M..]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  ans = {E(-2.29), ".5[[[[[[...]]]3mn[[....]]]MP]]"};
  EXPECT_EQ(ans, Mfe(m, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));

  // Has multiple possible MFE structures.
  EXPECT_EQ(E(-8.53), std::get<0>(Mfe(m, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG")));

  ans = {E(-19.99), "[[[[[Mn[[[[[[...]]]]]]]mp[[[[[[.......]]]]]]]M.m]]]]]3..."};
  EXPECT_EQ(ans, Mfe(m, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
  ans = {E(-88.99),
      "......m[[...[[[[[[[[............[[..[[[[[[[M..m[[[[........]]]]M.......[[[[[[n[[[[...[[.[[[."
      "....[[[.[[.[[[[......]]]].]].]]]......]]].]]]]]]]mp[[.[[[[....]]]]..]]]M...n[[[[[[[[[[......"
      ".......]]]]]]]]]]]mp[[[.....]]]]M...mn[...[[[[[[[.......]]]]]]]]]MP]]]]]3...........n[[[[[.."
      ".......]]]]]]mP]]]]]]..]]]]]]]]]]...]]M"};
  EXPECT_EQ(ans, Mfe(m, std::get<Primary>(k16sHSapiens3)));
}

#endif

INSTANTIATE_TEST_SUITE_P(MfeTest, MfeTestT12,
    testing::Combine(
        testing::Range(0, NUM_T12_MODELS), testing::ValuesIn(EnumValues<CtxCfg::MfeAlg>())));

}  // namespace mrna
