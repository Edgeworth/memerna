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

namespace mrna::md::t22 {

class T22MfeTest : public testing::TestWithParam<CtxCfg::DpAlg> {
 public:
  static std::tuple<Energy, std::string> Mfe(const erg::EnergyModelPtr& em, const std::string& s) {
    return Mfe(em, Primary::FromSeq(s));
  }

  static std::tuple<Energy, std::string> Mfe(const erg::EnergyModelPtr& em, const Primary& r) {
    auto res = Ctx(em, CtxCfg{.dp_alg = GetParam()}).Fold(r, {});
    return {res.mfe.energy, res.tb.ctd.ToString(res.tb.s)};
  }

  static void TestMfePseudofree(
      const Model::Ptr& base_em, Energy base_energy, const std::string& r, const std::string& db) {
    std::vector<Energy> pf_paired(r.size(), E(1.0));
    std::vector<Energy> pf_unpaired(r.size(), E(1.0));
    Energy extra_from_pseudofree = E(1.0) * int(r.size());
    auto em = base_em->CloneWithPseudofreeEnergy(pf_paired, pf_unpaired);
    auto [mfe_energy, mfe_db] = Mfe(em, r);
    EXPECT_EQ(base_energy + extra_from_pseudofree, mfe_energy);
    EXPECT_EQ(db, mfe_db);
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

TEST_P(T22MfeTest, T22P2PseudofreeEnergy) {
  auto em = t22p2;

  // Fast enough for brute force:
  TestMfePseudofree(em, E(-0.58), "CCUCCGGG", "[[....]]");
  TestMfePseudofree(em, E(-0.53), "CGGAAACGG", "[[....]]3");
  TestMfePseudofree(em, E(-0.49), "UGCAAAGCAA", "[[[...]]]3");
  TestMfePseudofree(em, E(-4.44), "GGGGAAACCCC", "[[[[...]]]]");
  TestMfePseudofree(em, E(-1.42), "CUUAUAGUUAAGG", "[[[[....]]]]3");
  TestMfePseudofree(em, E(-4.02), "CCGAAGGGGCUGCGGCG", "n[....]]p[....]].");
  TestMfePseudofree(em, E(-2.94), "GCCAAGGCCCCACCCGGA", "n[....]]mp[....]]M");
  TestMfePseudofree(em, E(-5.00), "GGCCGAUGGCAGCGAUAGC", "mn[[...]]]Mp[....]]");
  TestMfePseudofree(em, E(-2.08), "CUGAAACUGGAAACAGAAAUG", "......[[[....]]]3....");

  // Too slow for brute force:
  if (GetParam() == CtxCfg::DpAlg::BRUTE) return;
  TestMfePseudofree(em, E(-5.31), "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA",
      ".mn[...[[[....]]]]]Mp[[.[[[[..[[................]]..]]]]]]].");
  TestMfePseudofree(em, E(-13.48),
      "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU",
      "....n[[[[...]]]]]p[[[[[............[[..[[[...]]]..]]............]]]]]]..");
  TestMfePseudofree(
      em, E(-5.83), "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC", "...m[[[[[[....]]]3mn[[.....]]]MP]]M");
  TestMfePseudofree(em, E(-12.39), "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG",
      "m[[[[[3n[....]]mp[[[[.......]]]]]M]]]]]M.");
  TestMfePseudofree(em, E(-7.54), "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG",
      ".....5[[[[........[[[[[......]]]]]..]]]]");
  TestMfePseudofree(
      em, E(-3.24), "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG", "[[[3...n[[...]]]mp[.....]]M]]]3");
  TestMfePseudofree(em, E(-12.06), "CCGGGCCAGCCCGCUCCUACGGGGGGUC", ".n[[[....]]]]p[[[....]]]]...");
  TestMfePseudofree(
      em, E(-7.57), "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG", "[[[[M..[[[......]]]3n[[....]]]mP]]]");
  TestMfePseudofree(
      em, E(-3.75), "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC", "[[[[.[[.[[[....]]]...]].].]]]3.");
  TestMfePseudofree(em, E(-6.47), "CGCAGGGUCGGACCCGGGAGAACCGCGA", "[[[3n[[[...]]]]p[.....]]]]]3");
  TestMfePseudofree(
      em, E(-6.72), "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA", ".[[[[..[[[[....]]]].........]]]]3");
  TestMfePseudofree(em, E(-11.73),
      "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC",
      "...[[[[[M...[[.[[.[[[[...]]]].]].]]3.n[[..[[......]]..]]]mP]]]]");
  TestMfePseudofree(em, E(-4.57), "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC",
      "[[[3..n[[[........]]]]p[[[[....]]]]]..........]]]3");
  TestMfePseudofree(
      em, E(-5.95), "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC", "..m[[[[[[.....]]]]]]M..[[....]]3");
  TestMfePseudofree(em, E(-30.43),
      "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
      "[[[[[[[[[[.[[[.[[.[[....[[.[[[[[[......]]].]]].]]]].]].]]]...]]]..]]]]]]]3...");
  TestMfePseudofree(
      em, E(-5.33), "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU", "......n[[[....]]]]p[[......]]]");
  TestMfePseudofree(em, E(-15.33), "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG",
      "[[[3....n[[[[[[[[[....]]]]...]]]]]]p[[....]]].]]]3");
  TestMfePseudofree(
      em, E(-4.24), "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC", "[[[..........]]]3..............");
  TestMfePseudofree(
      em, E(-3.59), "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG", "[[[3n[[[....]]]]mp[[....]]]M...]]]3");
  TestMfePseudofree(
      em, E(-3.00), "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG", ".5[[[..[[[[....]]]]........]]]");
  TestMfePseudofree(
      em, E(-8.44), "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG", "[[[[[3n[[....]]]p[[[...]]]]]]]]]");
  TestMfePseudofree(em, E(-19.97), "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU",
      "[[[[[Mn[[[[[[...]]]]]]]mp[[[[[[.......]]]]]]]M.m]]]]]3...");
  TestMfePseudofree(em, E(-92.19), std::get<Primary>(k16sHSapiens3).ToSeq(),
      "......m[[...[[[[[[[[............[[..[[[[[[[Mp[.[[[[........]]]]..]]....n[[[[[n[[[[...[[.[[[."
      "....[[[.[[.[[[[......]]]].]].]]]......]]].]]]]]]]mp[[.[[[[....]]]]..]]]M...n[[[[[[[[[[......"
      ".......]]]]]]]]]]]mp[[[.....]]]]M...mn[...[[[[[[[.......]]]]]]]]]MP]]]]]p[[.[[[[....[[[...]]"
      "]....]]]].]]]mN]]]]]]..]]]]]]]]]]...]]M");
}

#else

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(T22MfeTest);

#endif

INSTANTIATE_TEST_SUITE_P(FoldAlgTest, T22MfeTest, testing::ValuesIn(CtxCfg::DP_ALGS));

}  // namespace mrna::md::t22
