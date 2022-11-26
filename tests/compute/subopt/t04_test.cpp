// Copyright 2022 Eliot Courtney.
#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

#include "common_test.h"
#include "compute/energy/model.h"
#include "compute/subopt/subopt.h"
#include "compute/subopt/subopt_cfg.h"
#include "ctx/ctx.h"
#include "ctx/ctx_cfg.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::subopt {

class SuboptAlgTest : public testing::TestWithParam<ctx::CtxCfg::SuboptAlg> {
 public:
  static std::vector<SuboptResult> Subopt(
      const erg::EnergyModelPtr& em, const std::string& s, const std::vector<Energy>& energies) {
    int n = static_cast<int>(energies.size());
    auto res = ctx::Ctx(em, ctx::CtxCfg{.subopt_alg = GetParam()})
                   .SuboptimalIntoVector(Primary::FromSeq(s), SuboptCfg{.strucs = n});
    for (int i = 0; i < n; ++i) EXPECT_EQ(res[i].energy, energies[i]);
    return res;
  }
};

#if ENERGY_PRECISION == 1

TEST_P(SuboptAlgTest, T04_P1) {
  Subopt(t04_p1, "CCUCCGGG",
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
  Subopt(t04_p1, "CGGAAACGG",
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
  Subopt(t04_p1, "UGCAAAGCAA",
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
  Subopt(t04_p1, "GGGGAAACCCC",
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
  Subopt(t04_p1, "GGGGAAACCCC",
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
  Subopt(t04_p1, "CUUAUAGUUAAGG",
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
  Subopt(t04_p1, "CCGAAGGGGCUGCGGCG",
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
  Subopt(t04_p1, "GCCAAGGCCCCACCCGGA",
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
  Subopt(t04_p1, "GGCCGAUGGCAGCGAUAGC",
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
  Subopt(t04_p1, "CUGAAACUGGAAACAGAAAUG",
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
  Subopt(t04_p1, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA",
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
  Subopt(t04_p1, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU",
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
  Subopt(t04_p1, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC",
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
  Subopt(t04_p1, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG",
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
  Subopt(t04_p1, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG",
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
  Subopt(t04_p1, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG",
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
  Subopt(t04_p1, "CCGGGCCAGCCCGCUCCUACGGGGGGUC",
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
  Subopt(t04_p1, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG",
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
  Subopt(t04_p1, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC",
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
  Subopt(t04_p1, "CGCAGGGUCGGACCCGGGAGAACCGCGA",
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
  Subopt(t04_p1, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA",
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
  Subopt(t04_p1, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC",
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
  Subopt(t04_p1, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC",
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
  Subopt(t04_p1, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC",
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
  Subopt(t04_p1, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
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
  Subopt(t04_p1, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU",
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
  Subopt(t04_p1, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG",
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
  Subopt(t04_p1, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC",
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
  Subopt(t04_p1, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG",
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
  Subopt(t04_p1, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG",
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
  Subopt(t04_p1, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG",
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
  Subopt(t04_p1, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU",
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

#elif ENERGY_PRECISION == 2

// TODO(0): Add tests for 2 decimal places.

#endif

INSTANTIATE_TEST_SUITE_P(SuboptAlgTest, SuboptAlgTest, testing::ValuesIn(ctx::CtxCfg::SUBOPT_ALGS));

}  // namespace mrna::subopt
