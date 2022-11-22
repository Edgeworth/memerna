// Copyright 2022 Eliot Courtney.
#include "compute/partition/partition.h"

#include <iomanip>
#include <string>

#include "common_test.h"
#include "compute/energy/energy.h"
#include "ctx/ctx.h"
#include "ctx/ctx_cfg.h"
#include "gtest/gtest.h"

namespace mrna::part {

class PartAlgTest : public testing::TestWithParam<ctx::CtxCfg::PartAlg> {
 public:
  static PartResult Partition(const std::string& s, energy::EnergyModelPtr em) {
    return ctx::Ctx(std::move(em), ctx::CtxCfg{.part_alg = GetParam()})
        .Partition(Primary::FromSeq(s));
  }
};

TEST_P(PartAlgTest, T04) {
  EXPECT_TRUE(rel_eq(4.24816013829, Partition("CCUCCGGG", t04).part.q));
  EXPECT_TRUE(rel_eq(4.17979557042, Partition("CGGAAACGG", t04).part.q));
  EXPECT_TRUE(rel_eq(4.62750397465, Partition("UGCAAAGCAA", t04).part.q));
  EXPECT_TRUE(rel_eq(3122.66115009, Partition("GGGGAAACCCC", t04).part.q));
  EXPECT_TRUE(rel_eq(9.44718000805, Partition("CUUAUAGUUAAGG", t04).part.q));
  EXPECT_TRUE(rel_eq(930.051311983, Partition("CCGAAGGGGCUGCGGCG", t04).part.q));
  EXPECT_TRUE(rel_eq(196.739899801, Partition("GCCAAGGCCCCACCCGGA", t04).part.q));
  EXPECT_TRUE(rel_eq(3431.19128513, Partition("GGCCGAUGGCAGCGAUAGC", t04).part.q));
  EXPECT_TRUE(rel_eq(94.3008892348, Partition("CUGAAACUGGAAACAGAAAUG", t04).part.q));

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::PartAlg::BRUTE) return;
  EXPECT_TRUE(rel_eq(573963557.833, Partition("CCGGGCCAGCCCGCUCCUACGGGGGGUC", t04).part.q));
  EXPECT_TRUE(rel_eq(226979.219097, Partition("CGCAGGGUCGGACCCGGGAGAACCGCGA", t04).part.q));
  EXPECT_TRUE(rel_eq(812.38352375, Partition("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG", t04).part.q));
  EXPECT_TRUE(rel_eq(19520.3202391, Partition("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU", t04).part.q));
  EXPECT_TRUE(rel_eq(2003.03838946, Partition("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG", t04).part.q));
  EXPECT_TRUE(rel_eq(23058.3336091, Partition("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC", t04).part.q));
  EXPECT_TRUE(rel_eq(2060.71972602, Partition("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC", t04).part.q));
  EXPECT_TRUE(rel_eq(10731115.7024, Partition("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG", t04).part.q));
  EXPECT_TRUE(rel_eq(443426.842362, Partition("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC", t04).part.q));
  EXPECT_TRUE(rel_eq(349171.530355, Partition("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA", t04).part.q));
  EXPECT_TRUE(rel_eq(132689.019045, Partition("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC", t04).part.q));
  EXPECT_TRUE(rel_eq(803546.641993, Partition("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG", t04).part.q));
  EXPECT_TRUE(rel_eq(2520.4316174, Partition("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG", t04).part.q));
  EXPECT_TRUE(
      rel_eq(643578.928228, Partition("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG", t04).part.q));
  EXPECT_TRUE(
      rel_eq(1660031715.98, Partition("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG", t04).part.q));
  EXPECT_TRUE(rel_eq(
      17711.6389272, Partition("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC", t04).part.q));
  EXPECT_TRUE(rel_eq(
      413455851608, Partition("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG", t04).part.q));
  EXPECT_TRUE(rel_eq(2.57353505099e+16,
      Partition("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU", t04).part.q));
  EXPECT_TRUE(rel_eq(265269.030237,
      Partition("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA", t04).part.q));
  EXPECT_TRUE(rel_eq(19869500930.9,
      Partition("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC", t04).part.q));
  EXPECT_TRUE(rel_eq(355492874401,
      Partition("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU", t04)
          .part.q));
  EXPECT_TRUE(rel_eq(1.40477999326e+21,
      Partition(
          "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA", t04)
          .part.q));
}

INSTANTIATE_TEST_SUITE_P(PartAlgTest, PartAlgTest, testing::ValuesIn(ctx::CtxCfg::PART_ALGS));

}  // namespace mrna::part
