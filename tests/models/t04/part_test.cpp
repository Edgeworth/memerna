// Copyright 2022 E.
#include "api/part.h"

#include <string>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/energy/model.h"
#include "common_test.h"
#include "gtest/gtest.h"
#include "model/primary.h"
#include "util/float.h"

namespace mrna::md::t04 {

class PartAlgTest : public testing::TestWithParam<CtxCfg::PartAlg> {
 public:
  static part::PartResult Partition(const std::string& s, const erg::EnergyModelPtr& em) {
    return Ctx(em, CtxCfg{.part_alg = GetParam()}).Partition(Primary::FromSeq(s));
  }
};

#if ENERGY_PRECISION == 1

TEST_P(PartAlgTest, T04P1) {
  EXPECT_TRUE(rel_eq(4.24816013829, Partition("CCUCCGGG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(4.17979557042, Partition("CGGAAACGG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(4.62750397465, Partition("UGCAAAGCAA", t04p1).part.q));
  EXPECT_TRUE(rel_eq(3122.66115009, Partition("GGGGAAACCCC", t04p1).part.q));
  EXPECT_TRUE(rel_eq(9.44718000805, Partition("CUUAUAGUUAAGG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(930.051311983, Partition("CCGAAGGGGCUGCGGCG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(196.739899801, Partition("GCCAAGGCCCCACCCGGA", t04p1).part.q));
  EXPECT_TRUE(rel_eq(3431.19128513, Partition("GGCCGAUGGCAGCGAUAGC", t04p1).part.q));
  EXPECT_TRUE(rel_eq(94.3008892348, Partition("CUGAAACUGGAAACAGAAAUG", t04p1).part.q));

  // Too slow for brute force:
  if (GetParam() == CtxCfg::PartAlg::BRUTE) return;
  EXPECT_TRUE(rel_eq(573963557.833, Partition("CCGGGCCAGCCCGCUCCUACGGGGGGUC", t04p1).part.q));
  EXPECT_TRUE(rel_eq(226979.219097, Partition("CGCAGGGUCGGACCCGGGAGAACCGCGA", t04p1).part.q));
  EXPECT_TRUE(rel_eq(812.38352375, Partition("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(19520.3202391, Partition("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU", t04p1).part.q));
  EXPECT_TRUE(rel_eq(2003.03838946, Partition("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(23058.3336091, Partition("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC", t04p1).part.q));
  EXPECT_TRUE(rel_eq(2060.71972602, Partition("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC", t04p1).part.q));
  EXPECT_TRUE(rel_eq(10731115.7024, Partition("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(443426.842362, Partition("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC", t04p1).part.q));
  EXPECT_TRUE(rel_eq(349171.530355, Partition("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA", t04p1).part.q));
  EXPECT_TRUE(
      rel_eq(132689.019045, Partition("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC", t04p1).part.q));
  EXPECT_TRUE(
      rel_eq(803546.641993, Partition("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(2520.4316174, Partition("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG", t04p1).part.q));
  EXPECT_TRUE(
      rel_eq(643578.928228, Partition("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG", t04p1).part.q));
  EXPECT_TRUE(
      rel_eq(1660031715.98, Partition("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(17711.6389272,
      Partition("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC", t04p1).part.q));
  EXPECT_TRUE(rel_eq(
      413455851608, Partition("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG", t04p1).part.q));
  EXPECT_TRUE(rel_eq(2.57353505099e+16,
      Partition("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU", t04p1).part.q));
  EXPECT_TRUE(rel_eq(265269.030237,
      Partition("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA", t04p1).part.q));
  EXPECT_TRUE(rel_eq(19869500930.9,
      Partition("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC", t04p1).part.q));
  EXPECT_TRUE(rel_eq(355492874401,
      Partition("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU", t04p1)
          .part.q));
  EXPECT_TRUE(rel_eq(1.40477999326e+21,
      Partition(
          "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA", t04p1)
          .part.q));
}

#elif ENERGY_PRECISION == 2

TEST_P(PartAlgTest, T04P2) {
  EXPECT_TRUE(rel_eq(4.06319224606, Partition("CCUCCGGG", t04p2).part.q));
  EXPECT_TRUE(rel_eq(3.99326566300, Partition("CGGAAACGG", t04p2).part.q));
  EXPECT_TRUE(rel_eq(5.00636522103, Partition("UGCAAAGCAA", t04p2).part.q));
  EXPECT_TRUE(rel_eq(2662.62888038, Partition("GGGGAAACCCC", t04p2).part.q));
  EXPECT_TRUE(rel_eq(10.7704965055, Partition("CUUAUAGUUAAGG", t04p2).part.q));
  EXPECT_TRUE(rel_eq(862.912629278, Partition("CCGAAGGGGCUGCGGCG", t04p2).part.q));
  EXPECT_TRUE(rel_eq(187.240327086, Partition("GCCAAGGCCCCACCCGGA", t04p2).part.q));
  EXPECT_TRUE(rel_eq(3427.00722809, Partition("GGCCGAUGGCAGCGAUAGC", t04p2).part.q));
  EXPECT_TRUE(rel_eq(93.019535283, Partition("CUGAAACUGGAAACAGAAAUG", t04p2).part.q));

  // Too slow for brute force:
  if (GetParam() == CtxCfg::PartAlg::BRUTE) return;
  EXPECT_TRUE(rel_eq(493466463.46, Partition("CCGGGCCAGCCCGCUCCUACGGGGGGUC", t04p2).part.q));
  EXPECT_TRUE(rel_eq(195295.14577, Partition("CGCAGGGUCGGACCCGGGAGAACCGCGA", t04p2).part.q));
  EXPECT_TRUE(rel_eq(778.887013804, Partition("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG", t04p2).part.q));
  EXPECT_TRUE(rel_eq(20972.1314181, Partition("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU", t04p2).part.q));
  EXPECT_TRUE(rel_eq(1909.72100698, Partition("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG", t04p2).part.q));
  EXPECT_TRUE(rel_eq(23557.9099628, Partition("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC", t04p2).part.q));
  EXPECT_TRUE(rel_eq(2064.1442172, Partition("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC", t04p2).part.q));
  EXPECT_TRUE(rel_eq(12061583.6231, Partition("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG", t04p2).part.q));
  EXPECT_TRUE(rel_eq(472928.685895, Partition("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC", t04p2).part.q));
  EXPECT_TRUE(rel_eq(365289.11509, Partition("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA", t04p2).part.q));
  EXPECT_TRUE(rel_eq(131826.62873, Partition("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC", t04p2).part.q));
  EXPECT_TRUE(
      rel_eq(870099.928286, Partition("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG", t04p2).part.q));
  EXPECT_TRUE(
      rel_eq(3100.78279704, Partition("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG", t04p2).part.q));
  EXPECT_TRUE(
      rel_eq(615865.695815, Partition("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG", t04p2).part.q));
  EXPECT_TRUE(
      rel_eq(1607489347.26, Partition("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG", t04p2).part.q));
  EXPECT_TRUE(rel_eq(23894.2123244,
      Partition("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC", t04p2).part.q));
  EXPECT_TRUE(rel_eq(
      421630654129, Partition("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG", t04p2).part.q));
  EXPECT_TRUE(rel_eq(2.41928077402e+16,
      Partition("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU", t04p2).part.q));
  EXPECT_TRUE(rel_eq(328476.821043,
      Partition("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA", t04p2).part.q));
  EXPECT_TRUE(rel_eq(17510378182.8,
      Partition("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC", t04p2).part.q));
  EXPECT_TRUE(rel_eq(389940780119,
      Partition("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU", t04p2)
          .part.q));
  EXPECT_TRUE(rel_eq(1.13634433186e+21,
      Partition(
          "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA", t04p2)
          .part.q));
}

// NEWMODEL: Add tests here.

#endif

INSTANTIATE_TEST_SUITE_P(PartAlgTest, PartAlgTest, testing::ValuesIn(CtxCfg::PART_ALGS));

}  // namespace mrna::md::t04
