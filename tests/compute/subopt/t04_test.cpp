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

TEST_P(SuboptAlgTest, T04P1) {
  Subopt(t04p1, "CCUCCGGG",
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
  Subopt(t04p1, "CGGAAACGG",
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
  Subopt(t04p1, "UGCAAAGCAA",
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
  Subopt(t04p1, "GGGGAAACCCC",
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
  Subopt(t04p1, "GGGGAAACCCC",
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
  Subopt(t04p1, "CUUAUAGUUAAGG",
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
  Subopt(t04p1, "CCGAAGGGGCUGCGGCG",
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
  Subopt(t04p1, "GCCAAGGCCCCACCCGGA",
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
  Subopt(t04p1, "GGCCGAUGGCAGCGAUAGC",
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
  Subopt(t04p1, "CUGAAACUGGAAACAGAAAUG",
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
  Subopt(t04p1, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA",
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
  Subopt(t04p1, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU",
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
  Subopt(t04p1, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC",
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
  Subopt(t04p1, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG",
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
  Subopt(t04p1, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG",
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
  Subopt(t04p1, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG",
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
  Subopt(t04p1, "CCGGGCCAGCCCGCUCCUACGGGGGGUC",
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
  Subopt(t04p1, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG",
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
  Subopt(t04p1, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC",
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
  Subopt(t04p1, "CGCAGGGUCGGACCCGGGAGAACCGCGA",
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
  Subopt(t04p1, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA",
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
  Subopt(t04p1, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC",
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
  Subopt(t04p1, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC",
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
  Subopt(t04p1, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC",
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
  Subopt(t04p1, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
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
  Subopt(t04p1, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU",
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
  Subopt(t04p1, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG",
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
  Subopt(t04p1, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC",
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
  Subopt(t04p1, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG",
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
  Subopt(t04p1, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG",
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
  Subopt(t04p1, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG",
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
  Subopt(t04p1, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU",
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

TEST_P(SuboptAlgTest, T04P2) {
  Subopt(t04p2, "CCUCCGGG",
      {
          E(-0.56),
          E(0.00),
          E(0.84),
          E(1.20),
          E(1.40),
          E(2.14),
          E(2.40),
          E(2.70),
          E(3.40),
          E(3.44),
      });
  Subopt(t04p2, "CGGAAACGG",
      {
          E(-0.56),
          E(0.00),
          E(0.74),
          E(1.40),
          E(1.60),
          E(2.60),
          E(2.80),
          E(3.10),
          E(3.90),
      });
  Subopt(t04p2, "UGCAAAGCAA",
      {
          E(-0.48),
          E(0.00),
          E(0.28),
          E(0.32),
          E(0.48),
          E(1.98),
          E(1.98),
          E(2.14),
          E(2.90),
          E(2.94),
      });
  Subopt(t04p2, "GGGGAAACCCC",
      {
          E(-4.38),
          E(-4.22),
          E(-3.42),
          E(-2.42),
          E(-2.42),
          E(-2.32),
          E(-2.12),
          E(-1.92),
          E(-1.37),
          E(-1.37),
      });
  Subopt(t04p2, "GGGGAAACCCC",
      {
          E(-4.38),
          E(-4.22),
          E(-3.42),
          E(-2.42),
          E(-2.42),
          E(-2.32),
          E(-2.12),
          E(-1.92),
          E(-1.37),
          E(-1.37),
      });
  Subopt(t04p2, "CUUAUAGUUAAGG",
      {
          E(-1.29),
          E(0.00),
          E(0.01),
          E(0.74),
          E(1.54),
          E(1.74),
          E(1.87),
          E(2.04),
          E(2.44),
          E(2.47),
      });
  Subopt(t04p2, "CCGAAGGGGCUGCGGCG",
      {
          E(-3.94),
          E(-2.38),
          E(-2.32),
          E(-2.22),
          E(-2.18),
          E(-1.50),
          E(-1.48),
          E(-1.42),
          E(-1.36),
          E(-1.12),
      });
  Subopt(t04p2, "GCCAAGGCCCCACCCGGA",
      {
          E(-2.88),
          E(-2.04),
          E(-1.18),
          E(-0.92),
          E(-0.92),
          E(-0.78),
          E(-0.72),
          E(-0.66),
          E(-0.49),
          E(-0.38),
      });
  Subopt(t04p2, "GGCCGAUGGCAGCGAUAGC",
      {
          E(-4.90),
          E(-3.10),
          E(-2.98),
          E(-2.80),
          E(-2.68),
          E(-2.34),
          E(-1.60),
          E(-1.60),
          E(-1.40),
          E(-1.40),
      });
  Subopt(t04p2, "CUGAAACUGGAAACAGAAAUG",
      {
          E(-2.19),
          E(-2.19),
          E(-1.59),
          E(-1.09),
          E(0.00),
          E(0.11),
          E(0.44),
          E(0.64),
          E(1.21),
          E(1.21),
      });

  // Too slow for brute force:
  if (GetParam() == ctx::CtxCfg::SuboptAlg::BRUTE) return;
  Subopt(t04p2, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA",
      {
          E(-5.25),
          E(-4.93),
          E(-4.93),
          E(-4.85),
          E(-4.80),
          E(-4.65),
          E(-4.64),
          E(-4.55),
          E(-4.53),
          E(-4.53),
      });
  Subopt(t04p2, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU",
      {
          E(-13.47),
          E(-13.37),
          E(-13.36),
          E(-13.17),
          E(-13.09),
          E(-13.06),
          E(-13.04),
          E(-12.99),
          E(-12.98),
          E(-12.98),
      });
  Subopt(t04p2, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC",
      {
          E(-5.69),
          E(-5.63),
          E(-5.49),
          E(-5.43),
          E(-5.43),
          E(-5.23),
          E(-5.19),
          E(-4.99),
          E(-4.95),
          E(-4.89),
      });
  Subopt(t04p2, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG",
      {
          E(-12.08),
          E(-11.78),
          E(-11.58),
          E(-11.28),
          E(-11.18),
          E(-10.98),
          E(-10.98),
          E(-10.88),
          E(-10.68),
          E(-10.68),
      });
  Subopt(t04p2, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG",
      {
          E(-7.33),
          E(-7.23),
          E(-6.37),
          E(-6.37),
          E(-5.79),
          E(-5.75),
          E(-5.75),
          E(-5.69),
          E(-5.65),
          E(-5.65),
      });
  Subopt(t04p2, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG",
      {
          E(-3.14),
          E(-2.74),
          E(-2.59),
          E(-2.54),
          E(-2.54),
          E(-2.45),
          E(-2.38),
          E(-2.34),
          E(-2.29),
          E(-2.15),
      });
  Subopt(t04p2, "CCGGGCCAGCCCGCUCCUACGGGGGGUC",
      {
          E(-11.90),
          E(-10.96),
          E(-10.78),
          E(-10.78),
          E(-10.12),
          E(-10.08),
          E(-10.00),
          E(-9.78),
          E(-9.62),
          E(-9.50),
      });
  Subopt(t04p2, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG",
      {
          E(-7.45),
          E(-7.15),
          E(-7.15),
          E(-6.85),
          E(-6.45),
          E(-6.25),
          E(-6.15),
          E(-5.97),
          E(-5.95),
          E(-5.85),
      });
  Subopt(t04p2, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC",
      {
          E(-2.97),
          E(-2.74),
          E(-2.54),
          E(-2.48),
          E(-2.46),
          E(-2.44),
          E(-2.28),
          E(-2.27),
          E(-2.26),
          E(-2.24),
      });
  Subopt(t04p2, "CGCAGGGUCGGACCCGGGAGAACCGCGA",
      {
          E(-6.41),
          E(-6.04),
          E(-5.65),
          E(-5.42),
          E(-5.37),
          E(-5.33),
          E(-5.33),
          E(-5.31),
          E(-5.30),
          E(-5.17),
      });
  Subopt(t04p2, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA",
      {
          E(-6.07),
          E(-5.89),
          E(-5.87),
          E(-5.83),
          E(-5.82),
          E(-5.69),
          E(-5.63),
          E(-5.63),
          E(-5.62),
          E(-5.59),
      });
  Subopt(t04p2, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC",
      {
          E(-12.04),
          E(-12.04),
          E(-12.04),
          E(-12.04),
          E(-11.92),
          E(-11.92),
          E(-11.92),
          E(-11.92),
          E(-11.44),
          E(-11.44),
      });
  Subopt(t04p2, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC",
      {
          E(-4.10),
          E(-4.00),
          E(-3.90),
          E(-3.85),
          E(-3.75),
          E(-3.71),
          E(-3.65),
          E(-3.55),
          E(-3.51),
          E(-3.46),
      });
  Subopt(t04p2, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC",
      {
          E(-6.73),
          E(-6.73),
          E(-6.67),
          E(-6.33),
          E(-6.33),
          E(-6.27),
          E(-6.03),
          E(-6.03),
          E(-5.97),
          E(-5.93),
      });
  Subopt(t04p2, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
      {
          E(-27.66),
          E(-27.10),
          E(-27.10),
          E(-26.86),
          E(-26.78),
          E(-26.74),
          E(-26.74),
          E(-26.72),
          E(-26.72),
          E(-26.69),
      });
  Subopt(t04p2, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU",
      {
          E(-5.35),
          E(-4.97),
          E(-4.87),
          E(-4.07),
          E(-3.99),
          E(-3.77),
          E(-3.74),
          E(-3.73),
          E(-3.67),
          E(-3.64),
      });
  Subopt(t04p2, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG",
      {
          E(-15.70),
          E(-15.50),
          E(-14.84),
          E(-14.64),
          E(-14.34),
          E(-14.20),
          E(-14.18),
          E(-14.14),
          E(-14.00),
          E(-14.00),
      });
  Subopt(t04p2, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC",
      {
          E(-4.42),
          E(-4.42),
          E(-4.22),
          E(-4.22),
          E(-4.18),
          E(-4.02),
          E(-4.02),
          E(-4.00),
          E(-3.92),
          E(-3.91),
      });
  Subopt(t04p2, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG",
      {
          E(-3.06),
          E(-2.86),
          E(-2.86),
          E(-2.86),
          E(-2.78),
          E(-2.58),
          E(-2.58),
          E(-2.58),
          E(-2.57),
          E(-2.40),
      });
  Subopt(t04p2, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG",
      {
          E(-2.30),
          E(-2.23),
          E(-1.94),
          E(-1.94),
          E(-1.88),
          E(-1.82),
          E(-1.80),
          E(-1.76),
          E(-1.76),
          E(-1.74),
      });
  Subopt(t04p2, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG",
      {
          E(-8.12),
          E(-8.12),
          E(-7.92),
          E(-7.92),
          E(-7.89),
          E(-7.89),
          E(-7.78),
          E(-7.70),
          E(-7.70),
          E(-7.59),
      });
  Subopt(t04p2, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU",
      {
          E(-20.82),
          E(-20.82),
          E(-20.82),
          E(-20.76),
          E(-20.76),
          E(-20.76),
          E(-20.56),
          E(-20.56),
          E(-20.56),
          E(-20.42),
      });
}

// NEWMODEL: Add tests here.

#endif

INSTANTIATE_TEST_SUITE_P(SuboptAlgTest, SuboptAlgTest, testing::ValuesIn(ctx::CtxCfg::SUBOPT_ALGS));

}  // namespace mrna::subopt
