// Copyright 2024 Eliot Courtney.
#include "api/subopt/subopt.h"

#include <string>
#include <tuple>
#include <vector>

#include "api/ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"
#include "tests/init.h"
#include "tests/util.h"

namespace mrna {

class SuboptTestT04 : public testing::TestWithParam<std::tuple<int, CtxCfg::SuboptAlg>> {
 public:
  static std::vector<subopt::SuboptResult> Subopt(
      const BackendModelPtr& m, const std::string& s, const std::vector<Energy>& energies) {
    return CheckSubopt(m, std::get<1>(GetParam()), s, energies);
  }

  static std::vector<subopt::SuboptResult> Subopt(
      const BackendModelPtr& m, const Primary& r, const std::vector<Energy>& energies) {
    return CheckSubopt(m, std::get<1>(GetParam()), r, energies);
  }
};

#if ENERGY_PRECISION == 1

TEST_P(SuboptTestT04, T04P1) {
  auto [i, alg] = GetParam();
  auto m = t04_ms[i];
  if (!Contains(CtxCfg::SuboptAlgsForBackend(m), alg)) return;

  Subopt(m, "CCUCCGGG",
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
  Subopt(m, "CGGAAACGG",
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
  Subopt(m, "UGCAAAGCAA",
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
  Subopt(m, "GGGGAAACCCC",
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
  Subopt(m, "GGGGAAACCCC",
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
  Subopt(m, "CUUAUAGUUAAGG",
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
  Subopt(m, "CCGAAGGGGCUGCGGCG",
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
  Subopt(m, "GCCAAGGCCCCACCCGGA",
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
  Subopt(m, "GGCCGAUGGCAGCGAUAGC",
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
  Subopt(m, "CUGAAACUGGAAACAGAAAUG",
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
  if (alg == CtxCfg::SuboptAlg::BRUTE) return;
  Subopt(m, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA",
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
  Subopt(m, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU",
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
  Subopt(m, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC",
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
  Subopt(m, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG",
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
  Subopt(m, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG",
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
  Subopt(m, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG",
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
  Subopt(m, "CCGGGCCAGCCCGCUCCUACGGGGGGUC",
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
  Subopt(m, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG",
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
  Subopt(m, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC",
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
  Subopt(m, "CGCAGGGUCGGACCCGGGAGAACCGCGA",
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
  Subopt(m, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA",
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
  Subopt(m, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC",
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
  Subopt(m, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC",
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
  Subopt(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC",
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
  Subopt(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
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
  Subopt(m, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU",
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
  Subopt(m, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG",
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
  Subopt(m, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC",
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
  Subopt(m, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG",
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
  Subopt(m, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG",
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
  Subopt(m, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG",
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
  Subopt(m, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU",
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
  Subopt(m, std::get<Primary>(k16sHSapiens3),
      {
          E(-89.6),
          E(-89.6),
          E(-89.6),
          E(-89.6),
          E(-89.5),
          E(-89.5),
          E(-89.5),
          E(-89.5),
          E(-89.5),
          E(-89.5),
      });
}

#elif ENERGY_PRECISION == 2

TEST_P(SuboptTestT04, T04P2) {
  auto [i, alg] = GetParam();
  const auto& m = t04_ms[i];
  if (!Contains(CtxCfg::SuboptAlgsForBackend(m), alg)) return;

  Subopt(m, "CCUCCGGG",
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
          E(3.49),
      });
  Subopt(m, "CGGAAACGG",
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
  Subopt(m, "UGCAAAGCAA",
      {
          E(-0.43),
          E(0.00),
          E(0.28),
          E(0.37),
          E(0.48),
          E(1.98),
          E(1.98),
          E(2.19),
          E(2.90),
          E(2.99),
      });
  Subopt(m, "GGGGAAACCCC",
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
  Subopt(m, "GGGGAAACCCC",
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
  Subopt(m, "CUUAUAGUUAAGG",
      {
          E(-1.24),
          E(0.00),
          E(0.06),
          E(0.79),
          E(1.64),
          E(1.84),
          E(1.92),
          E(2.09),
          E(2.54),
          E(2.57),
      });
  Subopt(m, "CCGAAGGGGCUGCGGCG",
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
  Subopt(m, "GCCAAGGCCCCACCCGGA",
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
  Subopt(m, "GGCCGAUGGCAGCGAUAGC",
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
  Subopt(m, "CUGAAACUGGAAACAGAAAUG",
      {
          E(-2.19),
          E(-2.19),
          E(-1.59),
          E(-1.09),
          E(0.00),
          E(0.11),
          E(0.49),
          E(0.69),
          E(1.21),
          E(1.31),
      });

  // Too slow for brute force:
  if (alg == CtxCfg::SuboptAlg::BRUTE) return;
  Subopt(m, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA",
      {
          E(-5.05),
          E(-4.73),
          E(-4.73),
          E(-4.65),
          E(-4.45),
          E(-4.45),
          E(-4.35),
          E(-4.33),
          E(-4.33),
          E(-4.29),
      });
  Subopt(m, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU",
      {
          E(-13.32),
          E(-13.27),
          E(-13.21),
          E(-13.07),
          E(-12.94),
          E(-12.91),
          E(-12.89),
          E(-12.89),
          E(-12.83),
          E(-12.83),
      });
  Subopt(m, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC",
      {
          E(-5.59),
          E(-5.58),
          E(-5.39),
          E(-5.38),
          E(-5.38),
          E(-5.18),
          E(-5.09),
          E(-4.89),
          E(-4.85),
          E(-4.84),
      });
  Subopt(m, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG",
      {
          E(-12.03),
          E(-11.73),
          E(-11.53),
          E(-11.23),
          E(-11.13),
          E(-10.93),
          E(-10.93),
          E(-10.83),
          E(-10.63),
          E(-10.63),
      });
  Subopt(m, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG",
      {
          E(-7.28),
          E(-7.18),
          E(-6.32),
          E(-6.32),
          E(-5.74),
          E(-5.70),
          E(-5.70),
          E(-5.64),
          E(-5.60),
          E(-5.60),
      });
  Subopt(m, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG",
      {
          E(-3.09),
          E(-2.69),
          E(-2.59),
          E(-2.49),
          E(-2.49),
          E(-2.45),
          E(-2.33),
          E(-2.29),
          E(-2.29),
          E(-2.15),
      });
  Subopt(m, "CCGGGCCAGCCCGCUCCUACGGGGGGUC",
      {
          E(-11.90),
          E(-10.91),
          E(-10.78),
          E(-10.78),
          E(-10.12),
          E(-10.08),
          E(-10.00),
          E(-9.78),
          E(-9.62),
          E(-9.45),
      });
  Subopt(m, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG",
      {
          E(-7.35),
          E(-7.05),
          E(-7.05),
          E(-6.75),
          E(-6.35),
          E(-6.15),
          E(-6.05),
          E(-5.92),
          E(-5.85),
          E(-5.75),
      });
  Subopt(m, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC",
      {
          E(-2.87),
          E(-2.64),
          E(-2.44),
          E(-2.43),
          E(-2.41),
          E(-2.34),
          E(-2.23),
          E(-2.22),
          E(-2.21),
          E(-2.18),
      });
  Subopt(m, "CGCAGGGUCGGACCCGGGAGAACCGCGA",
      {
          E(-6.36),
          E(-6.04),
          E(-5.65),
          E(-5.42),
          E(-5.32),
          E(-5.30),
          E(-5.28),
          E(-5.28),
          E(-5.26),
          E(-5.12),
      });
  Subopt(m, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA",
      {
          E(-5.97),
          E(-5.79),
          E(-5.77),
          E(-5.77),
          E(-5.73),
          E(-5.59),
          E(-5.57),
          E(-5.53),
          E(-5.53),
          E(-5.49),
      });
  Subopt(m, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC",
      {
          E(-11.99),
          E(-11.99),
          E(-11.99),
          E(-11.99),
          E(-11.82),
          E(-11.82),
          E(-11.82),
          E(-11.82),
          E(-11.39),
          E(-11.39),
      });
  Subopt(m, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC",
      {
          E(-3.95),
          E(-3.80),
          E(-3.75),
          E(-3.70),
          E(-3.60),
          E(-3.56),
          E(-3.50),
          E(-3.40),
          E(-3.36),
          E(-3.31),
      });
  Subopt(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC",
      {
          E(-6.68),
          E(-6.68),
          E(-6.62),
          E(-6.28),
          E(-6.28),
          E(-6.22),
          E(-5.98),
          E(-5.98),
          E(-5.92),
          E(-5.88),
      });
  Subopt(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
      {
          E(-27.41),
          E(-26.95),
          E(-26.95),
          E(-26.61),
          E(-26.59),
          E(-26.59),
          E(-26.58),
          E(-26.57),
          E(-26.57),
          E(-26.44),
      });
  Subopt(m, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU",
      {
          E(-5.25),
          E(-4.92),
          E(-4.82),
          E(-4.02),
          E(-3.94),
          E(-3.72),
          E(-3.68),
          E(-3.64),
          E(-3.62),
          E(-3.61),
      });
  Subopt(m, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG",
      {
          E(-15.65),
          E(-15.45),
          E(-14.79),
          E(-14.59),
          E(-14.29),
          E(-14.15),
          E(-14.09),
          E(-14.08),
          E(-13.95),
          E(-13.95),
      });
  Subopt(m, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC",
      {
          E(-4.42),
          E(-4.42),
          E(-4.22),
          E(-4.22),
          E(-4.18),
          E(-4.02),
          E(-4.02),
          E(-4.00),
          E(-3.91),
          E(-3.87),
      });
  Subopt(m, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG",
      {
          E(-2.91),
          E(-2.71),
          E(-2.71),
          E(-2.71),
          E(-2.63),
          E(-2.43),
          E(-2.43),
          E(-2.43),
          E(-2.42),
          E(-2.30),
      });
  Subopt(m, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG",
      {
          E(-2.15),
          E(-2.13),
          E(-1.84),
          E(-1.84),
          E(-1.83),
          E(-1.82),
          E(-1.71),
          E(-1.71),
          E(-1.65),
          E(-1.64),
      });
  Subopt(m, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG",
      {
          E(-8.07),
          E(-8.07),
          E(-7.89),
          E(-7.87),
          E(-7.87),
          E(-7.84),
          E(-7.73),
          E(-7.60),
          E(-7.60),
          E(-7.59),
      });
  Subopt(m, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU",
      {
          E(-20.62),
          E(-20.62),
          E(-20.62),
          E(-20.61),
          E(-20.61),
          E(-20.61),
          E(-20.41),
          E(-20.41),
          E(-20.41),
          E(-20.22),
      });
  Subopt(m, std::get<Primary>(k16sHSapiens3),
      {
          E(-89.29),
          E(-89.29),
          E(-89.26),
          E(-89.26),
          E(-89.19),
          E(-89.19),
          E(-89.16),
          E(-89.16),
          E(-89.13),
          E(-89.13),
      });
}

#endif

INSTANTIATE_TEST_SUITE_P(SuboptTest, SuboptTestT04,
    testing::Combine(
        testing::Range(0, NUM_T04_MODELS), testing::ValuesIn(EnumValues<CtxCfg::SuboptAlg>())));

}  // namespace mrna
