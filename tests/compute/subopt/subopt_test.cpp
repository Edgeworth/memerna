// Copyright 2022 Eliot Courtney.
#include "compute/subopt/subopt.h"

#include <string>
#include <vector>

#include "common_test.h"
#include "compute/subopt/config.h"
#include "ctx/config.h"
#include "ctx/ctx.h"
#include "gtest/gtest.h"
#include "model/primary.h"

namespace mrna::subopt {

class SuboptAlgTest : public testing::TestWithParam<ctx::CtxCfg::SuboptAlg> {
 public:
  static std::vector<SuboptResult> Subopt(const std::string& s) {
    return ctx::Ctx(t04, ctx::CtxCfg{.subopt_alg = GetParam()})
        .SuboptimalIntoVector(Primary::FromSeq(s), SuboptCfg{.strucs = 10});
  }
};

TEST_P(SuboptAlgTest, T04) {
  // TODO: Test these results. For now just run them.
  Subopt("CCUCCGGG");
  Subopt("CGGAAACGG");
  Subopt("UGCAAAGCAA");
  Subopt("GGGGAAACCCC");
  Subopt("GGGGAAACCCC");
  Subopt("CUUAUAGUUAAGG");
  Subopt("CCGAAGGGGCUGCGGCG");
  Subopt("GCCAAGGCCCCACCCGGA");
  Subopt("GGCCGAUGGCAGCGAUAGC");
  Subopt("CUGAAACUGGAAACAGAAAUG");

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::SuboptAlg::BRUTE) return;
  Subopt("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA");
  Subopt("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU");
  Subopt("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC");
  Subopt("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG");
  Subopt("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG");
  Subopt("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG");
  Subopt("CCGGGCCAGCCCGCUCCUACGGGGGGUC");
  Subopt("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG");
  Subopt("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC");
  Subopt("CGCAGGGUCGGACCCGGGAGAACCGCGA");
  Subopt("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA");
  Subopt("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC");
  Subopt("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC");
  Subopt("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC");
  Subopt("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA");
  Subopt("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU");
  Subopt("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG");
  Subopt("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC");
  Subopt("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG");
  Subopt("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG");
  Subopt("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG");
  Subopt("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU");
}

INSTANTIATE_TEST_SUITE_P(SuboptAlgTest, SuboptAlgTest, testing::ValuesIn(ctx::CtxCfg::SUBOPT_ALGS));

}  // namespace mrna::subopt
