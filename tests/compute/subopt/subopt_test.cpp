// Copyright 2022 Eliot Courtney.
#include "compute/subopt/subopt.h"

#include <sys/types.h>

#include <string>
#include <vector>

#include "common_test.h"
#include "compute/energy/model.h"
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
    constexpr auto NUM_CHECK = 10;
    verify(energies.size() == NUM_CHECK, "wrong number of energies");
    auto res = ctx::Ctx(em, ctx::CtxCfg{.subopt_alg = GetParam()})
                   .SuboptimalIntoVector(Primary::FromSeq(s), SuboptCfg{.strucs = NUM_CHECK});
    for (int i = 0; i < NUM_CHECK; ++i) EXPECT_EQ(res[i].energy, energies[i]);
    return res;
  }
};

TEST_P(SuboptAlgTest, T04) {
  Subopt("CCUCCGGG", t04,
      {
          E(0.60),
          E(0.00),
          E(0.80),
          E(1.20),
          E(1.40),
          E(2.10),
          E(2.40),
          E(2.70),
          E(3.40),
          E(3.50),
      });
  Subopt("CGGAAACGG", t04,
      {
          E(0.60),
          E(0.00),
          E(0.70),
          E(1.40),
          E(1.60),
          E(2.60),
          E(2.80),
          E(3.10),
          E(3.90),
      });
  Subopt("UGCAAAGCAA", t04,
      {
          E(0.40),
          E(0.00),
          E(0.30),
          E(0.40),
          E(0.50),
          E(2.00),
          E(2.00),
          E(2.20),
          E(2.90),
          E(3.00),
      });
  Subopt("GGGGAAACCCC", t04,
      {
          E(-4.50),
          E(-4.30),
          E(-3.50),
          E(-2.50),
          E(-2.50),
          E(-2.40),
          E(-2.20),
          E(-2.00),
          E(-1.50),
          E(-1.50),
      });
  Subopt("GGGGAAACCCC", t04,
      {
          E(-4.50),
          E(-4.30),
          E(-3.50),
          E(-2.50),
          E(-2.50),
          E(-2.40),
          E(-2.20),
          E(-2.00),
          E(-1.50),
          E(-1.50),
      });
  Subopt("CUUAUAGUUAAGG", t04,
      {
          E(-1.20),
          E(0.00),
          E(0.10),
          E(0.80),
          E(1.70),
          E(1.90),
          E(1.90),
          E(2.10),
          E(2.60),
          E(2.60),
      });
  Subopt("CCGAAGGGGCUGCGGCG", t04,
      {
          E(-4.00),
          E(-2.40),
          E(-2.30),
          E(-2.20),
          E(-2.20),
          E(-1.60),
          E(-1.50),
          E(-1.40),
          E(-1.40),
          E(-1.20),
      });
  Subopt("GCCAAGGCCCCACCCGGA", t04,
      {
          E(-2.90),
          E(-2.10),
          E(-1.20),
          E(-1.00),
          E(0.90),
          E(0.80),
          E(0.80),
          E(0.70),
          E(0.60),
          E(0.40),
      });
  Subopt("GGCCGAUGGCAGCGAUAGC", t04,
      {
          E(-4.90),
          E(-3.10),
          E(-3.00),
          E(-2.80),
          E(-2.70),
          E(-2.30),
          E(-1.60),
          E(-1.60),
          E(-1.40),
          E(-1.40),
      });
  Subopt("CUGAAACUGGAAACAGAAAUG", t04,
      {
          E(-2.20),
          E(-2.20),
          E(-1.60),
          E(-1.10),
          E(0.00),
          E(0.10),
          E(0.50),
          E(0.70),
          E(1.20),
          E(1.30),
      });

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::SuboptAlg::BRUTE) return;
  Subopt("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA", t04,
      {
          E(-5.10),
          E(-4.80),
          E(-4.80),
          E(-4.70),
          E(-4.50),
          E(-4.50),
          E(-4.40),
          E(-4.40),
          E(-4.40),
          E(-4.30),
      });
  Subopt("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU", t04,
      {
          E(-13.30),
          E(-13.30),
          E(-13.20),
          E(-13.10),
          E(-13.00),
          E(-13.00),
          E(-12.90),
          E(-12.90),
          E(-12.90),
          E(-12.90),
      });
  Subopt("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC", t04,
      {
          E(-5.70),
          E(-5.60),
          E(-5.50),
          E(-5.40),
          E(-5.40),
          E(-5.20),
          E(-5.20),
          E(-5.00),
          E(-5.00),
          E(-4.90),
      });
  Subopt("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG", t04,
      {
          E(-12.10),
          E(-11.80),
          E(-11.60),
          E(-11.30),
          E(-11.20),
          E(-11.00),
          E(-11.00),
          E(-10.90),
          E(-10.70),
          E(-10.70),
      });
  Subopt("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG", t04,
      {
          E(-7.40),
          E(-7.30),
          E(-6.30),
          E(-6.30),
          E(-5.80),
          E(-5.80),
          E(-5.80),
          E(-5.70),
          E(-5.70),
          E(-5.70),
      });
  Subopt("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG", t04,
      {
          E(-3.20),
          E(-2.80),
          E(-2.60),
          E(-2.60),
          E(-2.60),
          E(-2.50),
          E(-2.40),
          E(-2.40),
          E(-2.30),
          E(-2.20),
      });
  Subopt("CCGGGCCAGCCCGCUCCUACGGGGGGUC", t04,
      {
          E(-12.00),
          E(-11.00),
          E(-10.90),
          E(-10.90),
          E(-10.20),
          E(-10.20),
          E(-10.10),
          E(-9.90),
          E(-9.70),
          E(-9.60),
      });
  Subopt("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG", t04,
      {
          E(-7.40),
          E(-7.10),
          E(-7.10),
          E(-6.80),
          E(-6.40),
          E(-6.20),
          E(-6.10),
          E(-5.90),
          E(-5.90),
          E(-5.80),
      });
  Subopt("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC", t04,
      {
          E(-3.00),
          E(-2.70),
          E(-2.50),
          E(-2.50),
          E(-2.40),
          E(-2.40),
          E(-2.30),
          E(-2.30),
          E(-2.30),
          E(-2.20),
      });
  Subopt("CGCAGGGUCGGACCCGGGAGAACCGCGA", t04,
      {
          E(-6.50),
          E(-6.20),
          E(-5.80),
          E(-5.60),
          E(-5.40),
          E(-5.40),
          E(-5.40),
          E(-5.40),
          E(-5.40),
          E(-5.20),
      });
  Subopt("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA", t04,
      {
          E(-6.00),
          E(-5.90),
          E(-5.80),
          E(-5.80),
          E(-5.80),
          E(-5.70),
          E(-5.60),
          E(-5.60),
          E(-5.60),
          E(-5.60),
      });
  Subopt("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC", t04,
      {
          E(-12.20),
          E(-12.20),
          E(-12.20),
          E(-12.20),
          E(-12.00),
          E(-12.00),
          E(-12.00),
          E(-12.00),
          E(-11.60),
          E(-11.60),
      });
  Subopt("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC", t04,
      {
          E(-3.90),
          E(-3.80),
          E(-3.70),
          E(-3.60),
          E(-3.60),
          E(-3.50),
          E(-3.40),
          E(-3.40),
          E(-3.30),
          E(-3.30),
      });
  Subopt("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC", t04,
      {
          E(-6.70),
          E(-6.70),
          E(-6.60),
          E(-6.30),
          E(-6.30),
          E(-6.20),
          E(-6.00),
          E(-6.00),
          E(-5.90),
          E(-5.90),
      });
  Subopt("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA", t04,
      {
          E(-27.60),
          E(-27.30),
          E(-27.30),
          E(-26.90),
          E(-26.90),
          E(-26.90),
          E(-26.90),
          E(-26.80),
          E(-26.80),
          E(-26.80),
      });
  Subopt("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU", t04,
      {
          E(-5.30),
          E(-4.90),
          E(-4.80),
          E(-4.00),
          E(-4.00),
          E(-3.70),
          E(-3.70),
          E(-3.70),
          E(-3.60),
          E(-3.60),
      });
  Subopt("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG", t04,
      {
          E(-15.70),
          E(-15.50),
          E(-14.80),
          E(-14.60),
          E(-14.30),
          E(-14.20),
          E(-14.20),
          E(-14.10),
          E(-14.00),
          E(-14.00),
      });
  Subopt("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC", t04,
      {
          E(-4.40),
          E(-4.40),
          E(-4.20),
          E(-4.20),
          E(-4.20),
          E(-4.10),
          E(-4.00),
          E(-4.00),
          E(-3.90),
          E(-3.90),
      });
  Subopt("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG", t04,
      {
          E(-2.90),
          E(-2.70),
          E(-2.70),
          E(-2.70),
          E(-2.60),
          E(-2.40),
          E(-2.40),
          E(-2.40),
          E(-2.40),
          E(-2.30),
      });
  Subopt("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG", t04,
      {
          E(-2.30),
          E(-2.20),
          E(-2.00),
          E(-2.00),
          E(-2.00),
          E(-1.90),
          E(-1.80),
          E(-1.80),
          E(-1.80),
          E(-1.80),
      });
  Subopt("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG", t04,
      {
          E(-8.00),
          E(-8.00),
          E(-7.90),
          E(-7.80),
          E(-7.80),
          E(-7.80),
          E(-7.70),
          E(-7.60),
          E(-7.50),
          E(-7.50),
      });
  Subopt("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU", t04,
      {
          E(-20.80),
          E(-20.80),
          E(-20.80),
          E(-20.80),
          E(-20.80),
          E(-20.80),
          E(-20.60),
          E(-20.60),
          E(-20.60),
          E(-20.40),
      });
}

INSTANTIATE_TEST_SUITE_P(SuboptAlgTest, SuboptAlgTest, testing::ValuesIn(ctx::CtxCfg::SUBOPT_ALGS));

}  // namespace mrna::subopt
