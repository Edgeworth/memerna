// Copyright 2022 E.
#include "compute/subopt/subopt.h"

#include <string>

#include "common_test.h"
#include "ctx/ctx.h"
#include "gtest/gtest.h"

namespace mrna::subopt {

class SuboptAlgTest : public testing::TestWithParam<ctx::CtxCfg::SuboptAlg> {
 public:
  std::vector<SuboptResult> Subopt(const std::string& s) {
    return ctx::Ctx(t04, ctx::CtxCfg{.subopt_alg = GetParam()})
        .SuboptimalIntoVector(Primary::FromString(s), SuboptCfg{.strucs = 10});
  }
};

TEST_P(SuboptAlgTest, T04) {
  // TODO: Test these results. For now just run them.
  Subopt("GGGGAAACCCC");
  Subopt("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA");
  Subopt("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU");
  Subopt("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC");
  Subopt("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG");
  Subopt("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG");
  Subopt("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG");
  Subopt("CCGAAGGGGCUGCGGCG");
  Subopt("CCGGGCCAGCCCGCUCCUACGGGGGGUC");
  Subopt("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG");
  Subopt("CCUCCGGG");
  Subopt("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC");
  Subopt("CGCAGGGUCGGACCCGGGAGAACCGCGA");
  Subopt("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA");
  Subopt("CGGAAACGG");
  Subopt("CUGAAACUGGAAACAGAAAUG");
  Subopt("CUUAUAGUUAAGG");
  Subopt("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC");
  Subopt("GCCAAGGCCCCACCCGGA");
  Subopt("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC");
  Subopt("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC");
  Subopt("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA");
  Subopt("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU");
  Subopt("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG");
  Subopt("GGCCGAUGGCAGCGAUAGC");
  Subopt("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC");
  Subopt("GGGGAAACCCC");
  Subopt("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG");
  Subopt("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG");
  Subopt("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG");
  Subopt("UGCAAAGCAA");
  Subopt("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU");
}

INSTANTIATE_TEST_SUITE_P(SuboptAlgTest, SuboptAlgTest, testing::ValuesIn(ctx::CtxCfg::SUBOPT_ALGS));

}  // namespace mrna::subopt
