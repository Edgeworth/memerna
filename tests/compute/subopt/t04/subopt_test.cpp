// Copyright 2022 Eliot Courtney.
#include "compute/subopt/subopt.h"

#include <sys/types.h>

#include <string>
#include <vector>

#include "common_test.h"
#include "compute/energy/energy.h"
#include "compute/subopt/subopt_cfg.h"
#include "ctx/ctx.h"
#include "ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/primary.h"

namespace mrna::subopt {

class SuboptAlgTest : public testing::TestWithParam<ctx::CtxCfg::SuboptAlg> {
 public:
  static std::vector<SuboptResult> Subopt(
      const std::string& s, const energy::EnergyModelPtr& em, const std::vector<Energy>& energies) {
    int n = static_cast<int>(energies.size());
    auto res = ctx::Ctx(em, ctx::CtxCfg{.subopt_alg = GetParam()})
                   .SuboptimalIntoVector(Primary::FromSeq(s), SuboptCfg{.strucs = n});
    for (int i = 0; i < n; ++i) EXPECT_EQ(res[i].energy, energies[i]);
    return res;
  }
};

TEST_P(SuboptAlgTest, T04) {
  Subopt("CCUCCGGG", t04,
      {
          E(-0.6),
          E(0.0),
          E(0.8),
          E(1.2),
          E(1.4),
          E(2.1),
          E(2.4),
          E(2.7),
          E(3.4),
          E(3.5),
      });
  Subopt("CGGAAACGG", t04,
      {
          E(-0.6),
          E(0.0),
          E(0.7),
          E(1.4),
          E(1.6),
          E(2.6),
          E(2.8),
          E(3.1),
          E(3.9),
      });
  Subopt("UGCAAAGCAA", t04,
      {
          E(-0.4),
          E(0.0),
          E(0.3),
          E(0.4),
          E(0.5),
          E(2.0),
          E(2.0),
          E(2.2),
          E(2.9),
          E(3.0),
      });
  Subopt("GGGGAAACCCC", t04,
      {
          E(-4.5),
          E(-4.3),
          E(-3.5),
          E(-2.5),
          E(-2.5),
          E(-2.4),
          E(-2.2),
          E(-2.0),
          E(-1.5),
          E(-1.5),
      });
  Subopt("GGGGAAACCCC", t04,
      {
          E(-4.5),
          E(-4.3),
          E(-3.5),
          E(-2.5),
          E(-2.5),
          E(-2.4),
          E(-2.2),
          E(-2.0),
          E(-1.5),
          E(-1.5),
      });
  Subopt("CUUAUAGUUAAGG", t04,
      {
          E(-1.2),
          E(0.0),
          E(0.1),
          E(0.8),
          E(1.7),
          E(1.9),
          E(1.9),
          E(2.1),
          E(2.6),
          E(2.6),
      });
  Subopt("CCGAAGGGGCUGCGGCG", t04,
      {
          E(-4.0),
          E(-2.4),
          E(-2.3),
          E(-2.2),
          E(-2.2),
          E(-1.6),
          E(-1.5),
          E(-1.4),
          E(-1.4),
          E(-1.2),
      });
  Subopt("GCCAAGGCCCCACCCGGA", t04,
      {
          E(-2.9),
          E(-2.1),
          E(-1.2),
          E(-1.0),
          E(-0.9),
          E(-0.8),
          E(-0.8),
          E(-0.7),
          E(-0.6),
          E(-0.4),
      });
  Subopt("GGCCGAUGGCAGCGAUAGC", t04,
      {
          E(-4.9),
          E(-3.1),
          E(-3.0),
          E(-2.8),
          E(-2.7),
          E(-2.3),
          E(-1.6),
          E(-1.6),
          E(-1.4),
          E(-1.4),
      });
  Subopt("CUGAAACUGGAAACAGAAAUG", t04,
      {
          E(-2.2),
          E(-2.2),
          E(-1.6),
          E(-1.1),
          E(0.0),
          E(0.1),
          E(0.5),
          E(0.7),
          E(1.2),
          E(1.3),
      });

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::SuboptAlg::BRUTE) return;
  Subopt("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA", t04,
      {
          E(-5.1),
          E(-4.8),
          E(-4.8),
          E(-4.7),
          E(-4.5),
          E(-4.5),
          E(-4.4),
          E(-4.4),
          E(-4.4),
          E(-4.3),
      });
  Subopt("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU", t04,
      {
          E(-13.3),
          E(-13.3),
          E(-13.2),
          E(-13.1),
          E(-13.0),
          E(-13.0),
          E(-12.9),
          E(-12.9),
          E(-12.9),
          E(-12.9),
      });
  Subopt("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC", t04,
      {
          E(-5.7),
          E(-5.6),
          E(-5.5),
          E(-5.4),
          E(-5.4),
          E(-5.2),
          E(-5.2),
          E(-5.0),
          E(-5.0),
          E(-4.9),
      });
  Subopt("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG", t04,
      {
          E(-12.1),
          E(-11.8),
          E(-11.6),
          E(-11.3),
          E(-11.2),
          E(-11.0),
          E(-11.0),
          E(-10.9),
          E(-10.7),
          E(-10.7),
      });
  Subopt("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG", t04,
      {
          E(-7.4),
          E(-7.3),
          E(-6.3),
          E(-6.3),
          E(-5.8),
          E(-5.8),
          E(-5.8),
          E(-5.7),
          E(-5.7),
          E(-5.7),
      });
  Subopt("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG", t04,
      {
          E(-3.2),
          E(-2.8),
          E(-2.6),
          E(-2.6),
          E(-2.6),
          E(-2.5),
          E(-2.4),
          E(-2.4),
          E(-2.3),
          E(-2.2),
      });
  Subopt("CCGGGCCAGCCCGCUCCUACGGGGGGUC", t04,
      {
          E(-12.0),
          E(-11.0),
          E(-10.9),
          E(-10.9),
          E(-10.2),
          E(-10.2),
          E(-10.1),
          E(-9.9),
          E(-9.7),
          E(-9.6),
      });
  Subopt("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG", t04,
      {
          E(-7.4),
          E(-7.1),
          E(-7.1),
          E(-6.8),
          E(-6.4),
          E(-6.2),
          E(-6.1),
          E(-5.9),
          E(-5.9),
          E(-5.8),
      });
  Subopt("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC", t04,
      {
          E(-3.0),
          E(-2.7),
          E(-2.5),
          E(-2.5),
          E(-2.4),
          E(-2.4),
          E(-2.3),
          E(-2.3),
          E(-2.3),
          E(-2.2),
      });
  Subopt("CGCAGGGUCGGACCCGGGAGAACCGCGA", t04,
      {
          E(-6.5),
          E(-6.2),
          E(-5.8),
          E(-5.6),
          E(-5.4),
          E(-5.4),
          E(-5.4),
          E(-5.4),
          E(-5.4),
          E(-5.2),
      });
  Subopt("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA", t04,
      {
          E(-6.0),
          E(-5.9),
          E(-5.8),
          E(-5.8),
          E(-5.8),
          E(-5.7),
          E(-5.6),
          E(-5.6),
          E(-5.6),
          E(-5.6),
      });
  Subopt("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC", t04,
      {
          E(-12.2),
          E(-12.2),
          E(-12.2),
          E(-12.2),
          E(-12.0),
          E(-12.0),
          E(-12.0),
          E(-12.0),
          E(-11.6),
          E(-11.6),
      });
  Subopt("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC", t04,
      {
          E(-3.9),
          E(-3.8),
          E(-3.7),
          E(-3.6),
          E(-3.6),
          E(-3.5),
          E(-3.4),
          E(-3.4),
          E(-3.3),
          E(-3.3),
      });
  Subopt("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC", t04,
      {
          E(-6.7),
          E(-6.7),
          E(-6.6),
          E(-6.3),
          E(-6.3),
          E(-6.2),
          E(-6.0),
          E(-6.0),
          E(-5.9),
          E(-5.9),
      });
  Subopt("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA", t04,
      {
          E(-27.6),
          E(-27.3),
          E(-27.3),
          E(-26.9),
          E(-26.9),
          E(-26.9),
          E(-26.9),
          E(-26.8),
          E(-26.8),
          E(-26.8),
      });
  Subopt("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU", t04,
      {
          E(-5.3),
          E(-4.9),
          E(-4.8),
          E(-4.0),
          E(-4.0),
          E(-3.7),
          E(-3.7),
          E(-3.7),
          E(-3.6),
          E(-3.6),
      });
  Subopt("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG", t04,
      {
          E(-15.7),
          E(-15.5),
          E(-14.8),
          E(-14.6),
          E(-14.3),
          E(-14.2),
          E(-14.2),
          E(-14.1),
          E(-14.0),
          E(-14.0),
      });
  Subopt("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC", t04,
      {
          E(-4.4),
          E(-4.4),
          E(-4.2),
          E(-4.2),
          E(-4.2),
          E(-4.1),
          E(-4.0),
          E(-4.0),
          E(-3.9),
          E(-3.9),
      });
  Subopt("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG", t04,
      {
          E(-2.9),
          E(-2.7),
          E(-2.7),
          E(-2.7),
          E(-2.6),
          E(-2.4),
          E(-2.4),
          E(-2.4),
          E(-2.4),
          E(-2.3),
      });
  Subopt("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG", t04,
      {
          E(-2.3),
          E(-2.2),
          E(-2.0),
          E(-2.0),
          E(-2.0),
          E(-1.9),
          E(-1.8),
          E(-1.8),
          E(-1.8),
          E(-1.8),
      });
  Subopt("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG", t04,
      {
          E(-8.0),
          E(-8.0),
          E(-7.9),
          E(-7.8),
          E(-7.8),
          E(-7.8),
          E(-7.7),
          E(-7.6),
          E(-7.5),
          E(-7.5),
      });
  Subopt("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU", t04,
      {
          E(-20.8),
          E(-20.8),
          E(-20.8),
          E(-20.8),
          E(-20.8),
          E(-20.8),
          E(-20.6),
          E(-20.6),
          E(-20.6),
          E(-20.4),
      });
}

INSTANTIATE_TEST_SUITE_P(SuboptAlgTest, SuboptAlgTest, testing::ValuesIn(ctx::CtxCfg::SUBOPT_ALGS));

}  // namespace mrna::subopt
