// Copyright 2016 Eliot Courtney.
#include <string>
#include <tuple>

#include "api/ctx/backend.h"
#include "api/ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"
#include "tests/init.h"
#include "tests/util.h"

namespace mrna {

class MfeTestT22 : public testing::TestWithParam<std::tuple<int, CtxCfg::MfeAlg>> {
 public:
  static std::tuple<Energy, std::string> Mfe(const BackendModelPtr& m, const std::string& s) {
    return GetMfe(m, std::get<1>(GetParam()), s);
  }

  static std::tuple<Energy, std::string> Mfe(const BackendModelPtr& m, const Primary& r) {
    return GetMfe(m, std::get<1>(GetParam()), r);
  }

  static void TestMfePseudofree(const BackendModelPtr& base_em, Energy base_energy,
      const std::string& r, const std::string& db) {
    std::vector<Energy> pf_paired(r.size(), E(1.0));
    std::vector<Energy> pf_unpaired(r.size(), E(1.0));
    Energy extra_from_pseudofree = E(1.0) * int(r.size());
    auto m = CloneBackend(base_em);
    LoadPseudofreeEnergy(m, pf_paired, pf_unpaired);
    auto [mfe_energy, mfe_db] = Mfe(m, r);
    EXPECT_EQ(base_energy + extra_from_pseudofree, mfe_energy);
    EXPECT_EQ(db, mfe_db);
  }
};

#if ENERGY_PRECISION == 2

TEST_P(MfeTestT22, T22P2) {
  auto [i, alg] = GetParam();
  auto m = t22_ms[i];
  if (!Contains(CtxCfg::MfeAlgsForBackend(m), alg)) return;

  // Fast enough for brute force:
  std::tuple<Energy, std::string> ans = {E(-0.58), "[[....]]"};
  EXPECT_EQ(ans, Mfe(m, "CCUCCGGG"));
  ans = {E(-0.53), "[[....]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGGAAACGG"));
  ans = {E(-0.49), "[[[...]]]3"};
  EXPECT_EQ(ans, Mfe(m, "UGCAAAGCAA"));
  ans = {E(-4.44), "[[[[...]]]]"};
  EXPECT_EQ(ans, Mfe(m, "GGGGAAACCCC"));
  ans = {E(-1.42), "[[[[....]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CUUAUAGUUAAGG"));
  ans = {E(-4.02), "n[....]]p[....]]."};
  EXPECT_EQ(ans, Mfe(m, "CCGAAGGGGCUGCGGCG"));
  ans = {E(-2.94), "n[....]]mp[....]]M"};
  EXPECT_EQ(ans, Mfe(m, "GCCAAGGCCCCACCCGGA"));
  ans = {E(-5.00), "mn[[...]]]Mp[....]]"};
  EXPECT_EQ(ans, Mfe(m, "GGCCGAUGGCAGCGAUAGC"));
  ans = {E(-2.08), "......[[[....]]]3...."};
  EXPECT_EQ(ans, Mfe(m, "CUGAAACUGGAAACAGAAAUG"));

  // Too slow for brute force:
  if (alg == CtxCfg::MfeAlg::BRUTE) return;
  ans = {E(-5.31), ".mn[...[[[....]]]]]Mp[[.[[[[..[[................]]..]]]]]]]."};
  EXPECT_EQ(ans, Mfe(m, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  ans = {E(-13.48), "....n[[[[...]]]]]p[[[[[............[[..[[[...]]]..]]............]]]]]].."};
  EXPECT_EQ(
      ans, Mfe(m, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  ans = {E(-5.83), "...m[[[[[[....]]]3mn[[.....]]]MP]]M"};
  EXPECT_EQ(ans, Mfe(m, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  ans = {E(-12.39), "m[[[[[3n[....]]mp[[[[.......]]]]]M]]]]]M."};
  EXPECT_EQ(ans, Mfe(m, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  ans = {E(-7.54), ".....5[[[[........[[[[[......]]]]]..]]]]"};
  EXPECT_EQ(ans, Mfe(m, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  ans = {E(-3.24), "[[[3...n[[...]]]mp[.....]]M]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  ans = {E(-12.06), ".n[[[....]]]]p[[[....]]]]..."};
  EXPECT_EQ(ans, Mfe(m, "CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  ans = {E(-7.57), "[[[[M..[[[......]]]3n[[....]]]mP]]]"};
  EXPECT_EQ(ans, Mfe(m, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  ans = {E(-3.75), "[[[[.[[.[[[....]]]...]].].]]]3."};
  EXPECT_EQ(ans, Mfe(m, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  ans = {E(-6.47), "[[[3n[[[...]]]]p[.....]]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  ans = {E(-6.72), ".[[[[..[[[[....]]]].........]]]]3"};
  EXPECT_EQ(ans, Mfe(m, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  ans = {E(-11.73), "...[[[[[M...[[.[[.[[[[...]]]].]].]]3.n[[..[[......]]..]]]mP]]]]"};
  EXPECT_EQ(ans, Mfe(m, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  ans = {E(-4.57), "[[[3..n[[[........]]]]p[[[[....]]]]]..........]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  ans = {E(-5.95), "..m[[[[[[.....]]]]]]M..[[....]]3"};
  EXPECT_EQ(ans, Mfe(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  ans = {
      E(-30.43), "[[[[[[[[[[.[[[.[[.[[....[[.[[[[[[......]]].]]].]]]].]].]]]...]]]..]]]]]]]3..."};
  EXPECT_EQ(
      ans, Mfe(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  ans = {E(-5.33), "......n[[[....]]]]p[[......]]]"};
  EXPECT_EQ(ans, Mfe(m, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  ans = {E(-15.33), "[[[3....n[[[[[[[[[....]]]]...]]]]]]p[[....]]].]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  ans = {E(-4.24), "[[[..........]]]3.............."};
  EXPECT_EQ(ans, Mfe(m, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  ans = {E(-3.59), "[[[3n[[[....]]]]mp[[....]]]M...]]]3"};
  EXPECT_EQ(ans, Mfe(m, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  ans = {E(-3.00), ".5[[[..[[[[....]]]]........]]]"};
  EXPECT_EQ(ans, Mfe(m, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  ans = {E(-8.44), "[[[[[3n[[....]]]p[[[...]]]]]]]]]"};
  EXPECT_EQ(ans, Mfe(m, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  ans = {E(-19.97), "[[[[[Mn[[[[[[...]]]]]]]mp[[[[[[.......]]]]]]]M.m]]]]]3..."};
  EXPECT_EQ(ans, Mfe(m, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
  ans = {E(-92.19),
      "......m[[...[[[[[[[[............[[..[[[[[[[Mp[.[[[[........]]]]..]]....n[[[[[n[[[[...[[.[[[."
      "....[[[.[[.[[[[......]]]].]].]]]......]]].]]]]]]]mp[[.[[[[....]]]]..]]]M...n[[[[[[[[[[......"
      ".......]]]]]]]]]]]mp[[[.....]]]]M...mn[...[[[[[[[.......]]]]]]]]]MP]]]]]p[[.[[[[....[[[...]]"
      "]....]]]].]]]mN]]]]]]..]]]]]]]]]]...]]M"};
  EXPECT_EQ(ans, Mfe(m, std::get<Primary>(k16sHSapiens3)));
}

TEST_P(MfeTestT22, T22P2PseudofreeEnergy) {
  auto [i, alg] = GetParam();
  auto m = t22_ms[i];
  if (!Contains(CtxCfg::MfeAlgsForBackend(m), alg)) return;

  // Fast enough for brute force:
  TestMfePseudofree(m, E(-0.58), "CCUCCGGG", "[[....]]");
  TestMfePseudofree(m, E(-0.53), "CGGAAACGG", "[[....]]3");
  TestMfePseudofree(m, E(-0.49), "UGCAAAGCAA", "[[[...]]]3");
  TestMfePseudofree(m, E(-4.44), "GGGGAAACCCC", "[[[[...]]]]");
  TestMfePseudofree(m, E(-1.42), "CUUAUAGUUAAGG", "[[[[....]]]]3");
  TestMfePseudofree(m, E(-4.02), "CCGAAGGGGCUGCGGCG", "n[....]]p[....]].");
  TestMfePseudofree(m, E(-2.94), "GCCAAGGCCCCACCCGGA", "n[....]]mp[....]]M");
  TestMfePseudofree(m, E(-5.00), "GGCCGAUGGCAGCGAUAGC", "mn[[...]]]Mp[....]]");
  TestMfePseudofree(m, E(-2.08), "CUGAAACUGGAAACAGAAAUG", "......[[[....]]]3....");

  // Too slow for brute force:
  if (alg == CtxCfg::MfeAlg::BRUTE) return;
  TestMfePseudofree(m, E(-5.31), "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA",
      ".mn[...[[[....]]]]]Mp[[.[[[[..[[................]]..]]]]]]].");
  TestMfePseudofree(m, E(-13.48),
      "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU",
      "....n[[[[...]]]]]p[[[[[............[[..[[[...]]]..]]............]]]]]]..");
  TestMfePseudofree(
      m, E(-5.83), "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC", "...m[[[[[[....]]]3mn[[.....]]]MP]]M");
  TestMfePseudofree(m, E(-12.39), "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG",
      "m[[[[[3n[....]]mp[[[[.......]]]]]M]]]]]M.");
  TestMfePseudofree(m, E(-7.54), "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG",
      ".....5[[[[........[[[[[......]]]]]..]]]]");
  TestMfePseudofree(
      m, E(-3.24), "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG", "[[[3...n[[...]]]mp[.....]]M]]]3");
  TestMfePseudofree(m, E(-12.06), "CCGGGCCAGCCCGCUCCUACGGGGGGUC", ".n[[[....]]]]p[[[....]]]]...");
  TestMfePseudofree(
      m, E(-7.57), "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG", "[[[[M..[[[......]]]3n[[....]]]mP]]]");
  TestMfePseudofree(
      m, E(-3.75), "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC", "[[[[.[[.[[[....]]]...]].].]]]3.");
  TestMfePseudofree(m, E(-6.47), "CGCAGGGUCGGACCCGGGAGAACCGCGA", "[[[3n[[[...]]]]p[.....]]]]]3");
  TestMfePseudofree(
      m, E(-6.72), "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA", ".[[[[..[[[[....]]]].........]]]]3");
  TestMfePseudofree(m, E(-11.73), "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC",
      "...[[[[[M...[[.[[.[[[[...]]]].]].]]3.n[[..[[......]]..]]]mP]]]]");
  TestMfePseudofree(m, E(-4.57), "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC",
      "[[[3..n[[[........]]]]p[[[[....]]]]]..........]]]3");
  TestMfePseudofree(
      m, E(-5.95), "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC", "..m[[[[[[.....]]]]]]M..[[....]]3");
  TestMfePseudofree(m, E(-30.43),
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
      "[[[[[[[[[[.[[[.[[.[[....[[.[[[[[[......]]].]]].]]]].]].]]]...]]]..]]]]]]]3...");
  TestMfePseudofree(
      m, E(-5.33), "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU", "......n[[[....]]]]p[[......]]]");
  TestMfePseudofree(m, E(-15.33), "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG",
      "[[[3....n[[[[[[[[[....]]]]...]]]]]]p[[....]]].]]]3");
  TestMfePseudofree(
      m, E(-4.24), "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC", "[[[..........]]]3..............");
  TestMfePseudofree(
      m, E(-3.59), "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG", "[[[3n[[[....]]]]mp[[....]]]M...]]]3");
  TestMfePseudofree(
      m, E(-3.00), "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG", ".5[[[..[[[[....]]]]........]]]");
  TestMfePseudofree(
      m, E(-8.44), "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG", "[[[[[3n[[....]]]p[[[...]]]]]]]]]");
  TestMfePseudofree(m, E(-19.97), "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU",
      "[[[[[Mn[[[[[[...]]]]]]]mp[[[[[[.......]]]]]]]M.m]]]]]3...");
  TestMfePseudofree(m, E(-92.19), std::get<Primary>(k16sHSapiens3).ToSeq(),
      "......m[[...[[[[[[[[............[[..[[[[[[[Mp[.[[[[........]]]]..]]....n[[[[[n[[[[...[[.[[[."
      "....[[[.[[.[[[[......]]]].]].]]]......]]].]]]]]]]mp[[.[[[[....]]]]..]]]M...n[[[[[[[[[[......"
      ".......]]]]]]]]]]]mp[[[.....]]]]M...mn[...[[[[[[[.......]]]]]]]]]MP]]]]]p[[.[[[[....[[[...]]"
      "]....]]]].]]]mN]]]]]]..]]]]]]]]]]...]]M");
}

#else

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(MfeTestT22);

#endif

INSTANTIATE_TEST_SUITE_P(MfeTest, MfeTestT22,
    testing::Combine(
        testing::Range(0, NUM_T22_MODELS), testing::ValuesIn(EnumValues<CtxCfg::MfeAlg>())));

}  // namespace mrna
